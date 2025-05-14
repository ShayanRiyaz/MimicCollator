'''
MimicCollator.py

Utilities for collating MIMIC III and IV waveform records into matfile sets.

Usage:
    from mimic_collator import MimicCollator
    collator = MimicCollator(config)
    collator.scan_waveform_directory(...)

Author: Shayan Riyaz
License: MIT

'''

from rich import print
import os,time,json, requests, wfdb
from scipy.io import savemat
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterator, Optional, Tuple, Dict, Any
from datetime import datetime, timezone
from urllib.request import urlretrieve
from get_typecast import get_typecast
import h5py

session = requests.Session()

def concat_records(group,record_cols):
    data = group[record_cols].values.flatten()
    kept = [str(x) for x in data if pd.notna(x) and str(x).strip()!='']
    return ";".join(kept)

def write_subject_record(h5_file, subj_id, data):
    grp = h5_file.require_group(f"subjects/{subj_id}")

    for section, section_data in data.items():
        subgrp = grp.require_group(section)

        for k, v in section_data.items():
            if isinstance(v, str):
                # Store strings as UTF-8 encoded fixed-length bytes
                dt = h5py.string_dtype(encoding='utf-8')
                subgrp.create_dataset(k, data=v, dtype=dt)

            elif isinstance(v, (int, float, np.integer, np.floating)):
                # Use dtype for platform-independent storage
                subgrp.create_dataset(k, data=v)

            elif isinstance(v, np.ndarray):
                # Compress + chunk for large waveform arrays
                subgrp.create_dataset(k, data=v, compression="gzip", chunks=True)

            else:
                raise TypeError(f"Unsupported type for {section}/{k}: {type(v)}")

class MimicCollator():
    def __init__(self, config,verbose = False,search_cache=True):
        self.bVerbose = verbose
        self.search_cache = []

        # Meta Data
        self.version_num = get_typecast(config, ["version_num"],cast_fn=int, default=0)
        self.ethnicity_extract = config.get("ethnicity_extract", False)
        self.mimic_num = get_typecast(config, ["mimic_info","mimic_num"], cast_fn=str, default="4",verbose=self.bVerbose)        
        self.custom_label = config.get("custom_label", -1) 
        self.custom_records = get_typecast(config, ["custom_records"], cast_fn=lambda v: pd.DataFrame(v), default=pd.DataFrame())
        self.ethnicity_data = get_typecast(config, ["ethnicity_data"], cast_fn= bool, default=False)

        # Subject/sample parameters
        self.num_subjects = self.remaining_subjects = get_typecast(config, ["num_subjects"], cast_fn=int, default=50)

        # File paths
        self.root_dir  = get_typecast(config, ["paths","local","root_folder"],cast_fn=Path, default=os.getcwd())
        self.root_dir.mkdir(parents=True, exist_ok=True)
        self.downloaded_files_dir = self.root_dir / "downloaded_files"
        self.out_dir = self.root_dir / f'dataset_v{self.version_num}' / f'subject_{self.num_subjects}_v{self.version_num}'
        self.out_dir.mkdir(parents=True,exist_ok=True)
        self.mimic_path  = get_typecast(config, ["mimic_info","mimic_path"],cast_fn=str, default="mimic4wdb/0.1.0/")
    

        ## Constants
        self.MAX_SET_SIZE = get_typecast(config,[ "max_set_size"], cast_fn=int, default=50)
        self.MIN_FREQ = get_typecast(config, ["min_freq"], cast_fn=int, default=60)
        self.START_MINUTES = 5
        self.END_MINUTES = get_typecast(config, ["min_minutes"], cast_fn=int, default=60)
        self.START_SAMPLES = int(60*self.START_MINUTES*self.MIN_FREQ)
        self.END_SAMPLES = int(( self.END_MINUTES*60*self.MIN_FREQ)+ self.START_SAMPLES)
           
        # Signal Labels
        self.ecg_labels = [sig.lower() for sig in get_typecast(config, ["required_signals","ecg_labels"],cast_fn=list(), default=['ii'])]
        self.ppg_labels = [sig.lower() for sig in get_typecast(config, ["required_signals","ppg_labels"],cast_fn=list(), default=['pleth'])]
        self.abp_labels = [sig.lower() for sig in  get_typecast(config, ["required_signals","ppg_labels"],cast_fn=list(), default=['abp','art'])]

        # Custom MIMIC-4 settings
        if self.mimic_num == "4":
            self.mimic_matched_path = 'https://physionet.org/files/mimic-iv-ecg/1.0/'
        
        self.signal_presence_dicts = []  # One dict per subject
        self.all_signals_seen = set()    # Expanding set of signal labels

    def collate_dataset(self, load_waveforms=True):
        ''' 
        -> COLLATE_DATASET
            - Inputs:

            - Ouputs:
        '''
        print(f"\n --- Downloading and collating MIMIC {self.mimic_num} matched waveform subset ---")
        self.setup_paths()

        # Prepare Available Records
        if self.mimic_num == "3":
            match_records = self.prepare_mimic3_record_list()
        elif self.mimic_num == "4":
            match_records = self.prepare_mimic4_record_list()

        # Update available record length
        if self.num_subjects > match_records.shape[0]:
            if self.bVerbose: print(f'Not enough subjects available: Setting num_subjects from {self.num_subjects} to {match_records.shape[0]}')
            self.num_subjects = match_records.shape[0]

         # Check for existing records mathching length and metadata_requirement
        if self.search_cache:
            csv_path, existing_n  = self.find_matching_records_info_file(self.num_subjects)
        else:
            csv_path = []
            existing_n = None
        if (existing_n is not None) and (self.num_subjects <= existing_n) and csv_path:
            df = pd.read_csv(csv_path)
            if existing_n > self.num_subjects:
                df = df.sample(n=self.num_subjects, random_state=42).reset_index(drop=True)
                if self.bVerbose: print(f"Loaded {existing_n}→sampled to {self.num_subjects} from {csv_path.name}")
            else: 
                if self.bVerbose: print(f"Loaded {existing_n} from {csv_path.name}")
        else:
            if self.bVerbose: print("No best‐records CSV found; extracting new mapping")
            df = self.extract_matching_record_info(match_records)
            out = Path(self.root_dir) / f"mimic{self.mimic_num}_{df.shape[0]}_best_records.csv"
            df.to_csv(out, index=False)
            if self.bVerbose: print(f"Wrote new {out.name}")

        if self.num_subjects > df.shape[0]:
            self.num_subjects = self.remaining_subjects = df.shape[0]
        if load_waveforms:
            self.scan_waveform_directory(match_records,record_info = df)
        return 


    def setup_paths(self):
        print("\n - Setting up parameters")
        if not os.path.exists(self.root_dir):
            raise FileNotFoundError(f"Root folder does not exist: {self.root_dir}")
        print("Working directory verified:", self.root_dir)





    def prepare_mimic3_record_list(self):
        df_categories = {}
        recs = pd.read_csv(f"https://physionet.org/files/{self.mimic_path}/RECORDS-waveforms",header=None, names=["path"])
        df_categories = pd.DataFrame({
            'matching_records':recs["path"].str.rpartition('/')[0] + '/',
            'subject_id': recs["path"].str.split("/", expand=True)[1].str.lstrip("p"),
            'filename': recs["path"].str.split("/", expand=True)[2]})
        df_categories = df_categories.drop_duplicates(subset='subject_id', keep='first')
        if not self.custom_records.empty: 
            ids_in_df2 = set(self.custom_records["subject_id"])
            df_categories['in_df2'] = df_categories['subject_id'].isin(ids_in_df2)
            any_matches = df_categories['in_df2'].any()
            df_categories['notes'] = str()
            if any_matches:
                label_map = self.custom_records.set_index('subject_id')['label'].to_dict()
                df_categories['notes'] = df_categories['subject_id'].map(label_map).fillna("").astype(str)
        df_categories.reset_index(inplace=True,drop=True)
        return df_categories
    




    def prepare_mimic4_record_list(self):
        df_categories = {}
        dtype_cols = {i: str for i in range(16,22)}
        machine_measurments_path = self.downloaded_files_dir / 'machine_measurements.csv'
        url      = f"{self.mimic_matched_path}/machine_measurements.csv"

        if not machine_measurments_path.exists():
            machine_measurments_path.parent.mkdir(parents=True, exist_ok=True)
            urlretrieve(url, machine_measurments_path)   # C-optimized download
        machine_measurements_db = pd.read_csv(machine_measurments_path, dtype=dtype_cols)
        recs = pd.read_csv(f"https://physionet.org/files/{self.mimic_path}/RECORDS",header=None, names=["path"])
    

        df_categories = pd.DataFrame({
            'matching_records':recs["path"],
            'subject_id':recs["path"].str.split("/", expand=True)[2].str.lstrip("p")
        })

        record_cols = [f"report_{i}" for i in range(18)]  
        def concat_records(group):
            data = group[record_cols].values.flatten()
            kept = [str(x) for x in data if pd.notna(x) and str(x).strip()!='']
            return ";".join(kept)
        
        ids_in_df2 = set(machine_measurements_db['subject_id'])
        df_categories['in_df2'] = df_categories['subject_id'].isin(ids_in_df2)
        any_matches = df_categories['in_df2'].any()
        df_categories["notes"] = ""
        if any_matches:
            mapping = (machine_measurements_db.groupby("subject_id").apply(concat_records).rename("all_records"))
            df_categories = df_categories[df_categories["subject_id"].isin(mapping.index)].copy()
            df_categories["all_records"] = df_categories["subject_id"].map(mapping)
        return df_categories
    



    def find_matching_records_info_file(self, desired_n: int) -> Tuple[Optional[Path], int]:
        """
        Find a JSON metadata file with compatible labels and duration, then return:
        - the corresponding CSV path
        - and the number of subjects in that set

        Criteria:
        - smallest 'subjects_returned' ≥ desired_n (and meets other criteria)
        - if none, fallback to largest < desired_n
        - return (None, 0) if nothing matches
        """
        
        pattern = Path(self.root_dir) / f"mimic{self.mimic_num}_*_best_records_metadata.json"
        meta_files = sorted(pattern.parent.glob(pattern.name))

        supersets, subsets = [], []

        for meta_path in meta_files:
            try:
                with open(meta_path, 'r') as fp:
                    meta = json.load(fp)
                    
                subjects_available = meta.get("subjects_returned", 0)
                total_minutes = meta.get("end_minutes", 0)
                ppg_labels = meta.get("ppg_labels", [])
                ecg_labels = meta.get("ecg_labels", [])
                abp_labels = meta.get("abp_labels", [])

                is_match = (
                    subjects_available >= desired_n and
                    total_minutes >= self.END_MINUTES and
                    ppg_labels == self.ppg_labels and
                    ecg_labels == self.ecg_labels and
                    abp_labels == self.abp_labels
                )

                target_list = supersets if is_match else subsets
                target_list.append((subjects_available, meta_path))
                meta.close()
            except Exception as e:
                if self.bVerbose:
                    print(f"Skipping {meta_path} due to error: {e}")

        selected = (
            min(supersets, key=lambda x: x[0]) if supersets else
            max(subsets, key=lambda x: x[0]) if subsets else
            (None, 0)
        )

        if selected[0] is None:
            return None, 0

        # Get corresponding CSV path from metadata file name
        _, meta_file = selected
        csv_file = meta_file.with_name(meta_file.name.replace("_metadata.json", ".csv"))

        if not csv_file.exists():
            if self.bVerbose:
                print(f"Found valid metadata at {meta_file} but missing CSV: {csv_file}")
            return None, 0

        return csv_file, selected[0]

    def update_signal_presence(self, actual_signals):
        """
        Updates self.signal_presence_dicts and self.all_signals_seen
        based on the current subject's actual_signals list.

        Ensures:
        - All previous subjects get 0 for newly seen signals
        - Current subject gets a 1-hot dict over all seen signals
        """
        # Normalize current signal names
        signal_set = set(s.lower().strip() for s in actual_signals)

        # Find new signal labels not seen before
        new_signals = signal_set - self.all_signals_seen

        # Expand the global signal set
        self.all_signals_seen.update(new_signals)

        # Backfill 0s for the new signals in all previous records
        for d in self.signal_presence_dicts:
            for sig in new_signals:
                d[sig] = 0

        # Create the current subject's signal label dict
        signal_dict = {sig: (1 if sig in signal_set else 0) for sig in self.all_signals_seen}

        # Append to history
        self.signal_presence_dicts.append(signal_dict)

        return signal_dict

    def signal_requirements_check(self,filename,pn_dir,header_type):
        """Read the .hea header, stripping out any '~' segments."""
        try:
            hdr = wfdb.rdheader(filename, pn_dir=pn_dir)
        except Exception as e:
            return None,0,None
            
        if header_type == 'folder':
            if hasattr(hdr, "seg_name") and "~" in hdr.seg_name:
                bad = [i for i,n in enumerate(hdr.seg_name) if n == "~"]
                hdr.seg_name = [n for i,n in enumerate(hdr.seg_name) if i not in bad]
                hdr.seg_len  = np.delete(hdr.seg_len, bad)
            if not (hdr.sig_name):
                return hdr,1,None
        
        try:
            actual_signals = [s.lower() for s in hdr.sig_name]
        except:
            actual_signals = None
            return hdr,0,actual_signals
        
        if header_type == 'file':
            if (hdr.fs < self.MIN_FREQ):
                return hdr,0,actual_signals
            
        if self.abp_labels is not None:
            arterial_signals_present = [label.lower() for label in self.abp_labels if label  in actual_signals]
            if not arterial_signals_present: 
                return hdr,0,actual_signals

        if self.ecg_labels is not None:
            ecg_labels_present = all(label.lower() in actual_signals for label in self.ecg_labels)
            if not ecg_labels_present: 
                return hdr,0,actual_signals

        if self.ppg_labels is not None:
            ppg_labels_present = all(label.lower() in actual_signals for label in self.ppg_labels)
            if not ppg_labels_present: 
                return hdr,0,actual_signals

        return hdr,1,actual_signals




    def _extract_one_record(self, idx, curr_record, records_cache, df_categories):
        """
        The logic you had inside your for‐loop, reduced to a function that
        takes (idx, curr_record) and either returns an info‐dict or None.
        """
        rec_id = records_cache.get(curr_record)
        notes = df_categories.loc[df_categories["matching_records"] == curr_record,"notes"].squeeze() or ""
        if not rec_id:
            return None

        subject_and_rec_dir = os.path.join(self.mimic_path, curr_record)
        if self.mimic_num == "3":
            base_pn_dir     = subject_and_rec_dir
            initial_header  = rec_id[:-2] if rec_id.endswith("n") else rec_id
        else:
            base_pn_dir = os.path.join(subject_and_rec_dir, rec_id)
            initial_header  = rec_id

        if self.mimic_num == "3" and initial_header.endswith("n"):
            # Pending Work add vital type signs available
            return None

        try:
            hdr1, ok1, _ = self.signal_requirements_check(filename=initial_header,pn_dir=base_pn_dir,header_type="folder")
        except:
            return None
        if not ok1:
            return None

        segs = list(zip(hdr1.seg_name, hdr1.seg_len))
        segs.sort(key=lambda x: x[1], reverse=True)
        if len(segs) < 2 or segs[0][1] < int(self.END_SAMPLES):
            return None

        # second requirements check
        try:
            hea_signal_prop = segs[0][0].split(".")[0]
            _, ok2, actual_signals = self.signal_requirements_check(filename=hea_signal_prop,pn_dir=base_pn_dir,header_type="folder")
        except:
            return None
        if not ok2:
            signal_dict = self.update_signal_presence(actual_signals)
            # if self.bVerbose:       
            #     print(
            #     f"rej subject:  "
            #     f"idx: {str(idx).rjust(3)}, "
            #     # f"pn_dir: {base_pn_dir.ljust(35)}, "
            #     # f"rec_id: {rec_id.ljust(30)}, "
            #     # f"file_id: {segs[0][0].ljust(15)}, "
            #     # # f"signals: {str(actual_signals).ljust(40)}, "
            #     # f"notes: {notes.ljust(50)}",
            #     f"signal_dict: {signal_dict}"
            #     )
            
            # return {
            #     "idx":      idx,
            #     "pn_dir":   subject_and_rec_dir,
            #     "rec_id":   str(rec_id),
            #     "file_id":  segs[0][0],
            #     "seg_len":  segs[0][1],
            #     "max_freq": 0,
            #     "signals":  actual_signals,
            #     "notes":    notes,
            #     "signal_labels": signal_dict
            # }
            return None

        # finally loop segments until one passes file‐type check
        for seg_name, seg_len in segs:
            if seg_len < int(self.END_SAMPLES):
                break
            try:
                hdr3, ok3, actual_signals = self.signal_requirements_check(filename=seg_name,pn_dir=base_pn_dir,header_type="file")
            except:
                continue
            if not ok3:
                continue
            
            signal_dict = self.update_signal_presence(actual_signals)

            if self.bVerbose:
                print(
                f"add subject:  "
                f"idx: {str(idx).rjust(3)}, "
                f"pn_dir: {base_pn_dir.ljust(35)}, "
                f"rec_id: {rec_id.ljust(30)}, "
                f"file_id: {seg_name.ljust(15)}, "
                f"signals: {str(actual_signals).ljust(40)}, "
                f"notes: {notes.ljust(50)}",
                # f"signal_dict: {signal_dict}"
                )
            return {
                "idx":      idx,
                "pn_dir":   subject_and_rec_dir,
                "rec_id":   str(rec_id),
                "file_id":  seg_name,
                "seg_len":  hdr3.sig_len,
                "max_freq": hdr3.fs,
                "signals":  actual_signals,
                "notes":    notes,
                "signal_labels": signal_dict
            }
        return None

    def extract_matching_record_info(self, df_categories, max_workers=8):
        matching_records = df_categories['matching_records']
        if self.mimic_num == "3":
            records_cache = df_categories.set_index('matching_records').to_dict()['filename']            
        elif self.mimic_num == "4":
            records_cache = {}
            for rec in matching_records:
                url = f"https://physionet.org/files/{self.mimic_path}{rec}RECORDS"
                txt = session.get(url).text
                folder = txt.splitlines()[0].rstrip("/").rsplit("/", 1)[-1]
                records_cache[rec] = folder

        if isinstance(matching_records, dict):
            record_list = list(matching_records.values())
        else:
            record_list = list(matching_records)

        # build args for each record
        args = [(idx, rec, records_cache, df_categories) for idx, rec in enumerate(record_list)]
        infos = []
        bDone = False
        num_subjects_observed = 0
        try:
            with ThreadPoolExecutor(max_workers=max_workers) as exe:
                futures = {
                    exe.submit(self._extract_one_record, *arg): arg
                    for arg in args
                }
                for fut in as_completed(futures):
                    info = fut.result()
                    if info is not None:
                        infos.append(info)
                        num_subjects_observed +=1
                        if num_subjects_observed >= self.num_subjects:
                            bDone = True
                            break
                        # Exit the executor early by shutting it down explicitly
                if bDone:
                    exe.shutdown(wait=False, cancel_futures=True)
        except KeyboardInterrupt:
            print("Interrupted. Cancelling futures...")
            for f in futures:
                f.cancel()
            raise
        df = pd.DataFrame(infos)
        df = df.sort_values(by=['idx'], ascending=True)

        meta = {
            "generated_at":  datetime.now(timezone.utc).isoformat() + "Z",
            "timestamp_unix": time.time(),
            "start_minutes": self.START_MINUTES,
            "end_minutes":   self.END_MINUTES,
            "ppg_labels":    self.ppg_labels,
            "ecg_labels":    self.ecg_labels,
            "abp_labels":    self.abp_labels,
            "mimic_num":     self.mimic_num,
            "subjects_requested":  self.num_subjects,
            "subjects_returned":   num_subjects_observed,
        }
        meta_path = os.path.join(
            self.root_dir,
            f"mimic{self.mimic_num}_{df.shape[0]}_best_records_metadata.json"
        )
        with open(meta_path, "w") as mf:
            json.dump(meta, mf, indent=2)
        mf.close()
        if self.bVerbose:
            print(f"Saved metadata → {meta_path}")
        return df
    
    def get_records_iterator(self,df_categories: Dict[str, Any],
        record_info: Optional[pd.DataFrame]) -> Iterator[Tuple[int,str]]:
        """
        Yield (idx, record_path) in the correct order, respecting record_info if given.
        """
        matching = df_categories["matching_records"]
        if record_info is not None:
            matching = matching.iloc[record_info["idx"]]
        if hasattr(matching, "items"):
            yield from matching.items()
        else:
            yield from enumerate(matching)





    def compute_paths(self, curr_record: str, curr_idx: int, record_info: Optional[pd.DataFrame]) -> Tuple[str, str, str, str]:
        """
        Returns (pn_dir, hea_file, rec_id, file_id).
        Encapsulates all the logic for figuring out filesystem vs URL paths.
        """
        base = os.path.join(self.mimic_path, curr_record)
        if record_info is not None:
            mask = record_info["idx"] == curr_idx
            rec_id = str(record_info.loc[mask, "rec_id"].iat[0])
            file_id = str(record_info.loc[mask, "file_id"].iat[0])
            hea_file = f"{file_id}.hea"
            if self.mimic_num == "4":
                pn_dir = os.path.join(base, rec_id)
            else:
                pn_dir = base
        else:
            # the old “download then read RECORDS” fallback
            data_url = f"https://physionet.org/files/{self.mimic_path}{curr_record}"
            temp_folder = pd.read_csv(f"{data_url}/RECORDS", header=None).iat[0, 0].rsplit("/", 1)[-1]
            hea_file = f"{temp_folder}.hea"
            pn_dir = os.path.join(base, temp_folder) if self.mimic_num == "4" else base
            rec_id = temp_folder
            file_id = temp_folder
        return pn_dir, hea_file, rec_id, file_id

    def build_subject_data(self, record, subj_id: str, rec_id: str,file_id: str, notes: Any) -> Dict[str,Any]:
        """Compose the dict for this one subject."""
        sigs = [s.lower() for s in record.sig_name]
        abp_l = next(l for l in self.abp_labels if l in sigs)

        return {
        "fix": {
            "subj_id":     subj_id,
            "rec_id":      rec_id,
            "files":       file_id,
            "af_status":   self.custom_label,
            "subject_notes": notes
        },
        "ppg": {
            "v":      record.p_signal[:, 0],
            "fs":     record.fs,
            "method": "pleth from .hea/.dat",
            "label" : "pleth"
        },
        "ekg": {
            "v":      record.p_signal[:, 1],
            "fs":     record.fs,
            "method": "ECG from lead II",
            "label":  "ii"
        },
        "bp": {
            "v":      record.p_signal[:, 2],
            "fs":     record.fs,
            "method": f"{abp_l} from .hea/.dat",
            "label":  abp_l 
            }
        }

    def process_single_record(self, args):
        """
        Worker for one (idx, curr_record) tuple.
        Returns a dict of subject_data or None on failure / skip.
        """
        idx, curr_record, df_categories, record_info = args
        record_info = record_info.sort_values(by=['idx'], ascending=True)
        try:
            pn_dir, hea_file, rec_id, file_id = self.compute_paths(
                curr_record, idx, record_info)

            file_header,pass_check,sigs = self.signal_requirements_check(filename=hea_file.split('.')[0],pn_dir=pn_dir,header_type='file')
            seg_name = file_header.record_name

            if pass_check == 0:
                return data
            ppg_i = [i for i, sig in enumerate(sigs) if sig.lower() in self.ppg_labels]
            ecg_i = [i for i, sig in enumerate(sigs) if sig.lower() in self.ecg_labels]
            abp_l = next(l for l in self.abp_labels if l in sigs)
            abp_i = sigs.index(abp_l)
            channels = [ppg_i,ecg_i,abp_i]
            channels = [x for sub in channels for x in (sub if isinstance(sub, list) else [sub])]
            record = wfdb.rdrecord(seg_name, pn_dir=pn_dir, sampfrom=int(self.START_MINUTES*60*file_header.fs), sampto=int(self.END_MINUTES*60*file_header.fs), channels = channels, return_res=32)

            subj_id     = curr_record.split("/")[-2]
            subj_notes  = df_categories.loc[df_categories['subject_id'] == subj_id[1:]]['notes'].values
            data        = self.build_subject_data( record, subj_id, rec_id, file_id, subj_notes)
            return data

        except Exception as e:
            if self.bVerbose:
                print(f"Skipping {curr_record!r}: {e}")
            return None

    def scan_waveform_directory(self, df_categories: Dict[str,Any], record_info: Optional[pd.DataFrame], max_workers:int = 8) -> list:
        """
        Parallelized scan over get_records_iterator.
        """

        # subject_data      = []
        bDone = False
        set_num           = 0
        num_subjects_added = 0
        num_remaining_subjects = self.num_subjects
        arg_iter = [
            (idx, rec, df_categories, record_info)
            for idx, rec in self.get_records_iterator(df_categories, record_info)
        ]
        with ThreadPoolExecutor(max_workers=max_workers) as exe:
            futures = {exe.submit(self.process_single_record, args): args for args in arg_iter}        
            for fut in as_completed(futures):
                data = fut.result()
                if data is None:
                    continue

                # subject_data.append(data)
                num_subjects_added += 1
                num_remaining_subjects -= 1
                if (num_subjects_added >= self.MAX_SET_SIZE) or (num_remaining_subjects == 0):
                    set_num += 1
                    fname = f"mimic{self.mimic_num}_data_{self.num_subjects}_{set_num}.h5"
                    out_path = os.path.join(self.out_dir, fname)

                    with h5py.File(out_path, "w") as h5f:
                        write_subject_record(h5f, data['fix']['subj_id'], data)

                    if self.bVerbose: print(f"Saved {num_subjects_added} Subjects → {fname}")
                    num_subjects_added = 0
                if self.bVerbose: print(f"Remainings Subject ({num_remaining_subjects}/{self.num_subjects-set_num*self.MAX_SET_SIZE}) | {num_subjects_added} Subjects Added in Set #{set_num}")
                
                if num_subjects_added >= self.num_subjects:
                    bDone = True
                    break
            if bDone:
                exe.shutdown(wait=False, cancel_futures=True)
        return #subject_data
    
    # In Development
    # def download_and_extract_gz(self, url, out_path):
    #     print(f"Downloading from {url}")
    #     gz_path = out_path + ".gz"
    #     with requests.get(url, stream=True) as r:
    #         r.raise_for_status()
    #         with open(gz_path, 'wb') as f:
    #             shutil.copyfileobj(r.raw, f)
    #     with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
    #         shutil.copyfileobj(f_in, f_out)
    #     os.remove(gz_path)

# if __name__ == "__main__":
    
#     VERSION_NUM = 1
#     MIN_MINUTES = 30
#     ecg_labels = ['II']
#     ppg_labels = ['Pleth']
#     abp_labels = ['abp','art']
#     custom_label = 0

#     for MIMIC_NUM in ["3","4"]:
#         root_folder = f'data/raw/mimic{MIMIC_NUM}_data/'
#         if not os.path.exists(root_folder):
#             os.makedirs(root_folder, exist_ok = False)

#         custom_records = {}
#         bCustomRecords = True
#         if MIMIC_NUM == "3":
#             MIMIC_PATH = 'mimic3wdb-matched/1.0/'
#             MAX_SUBJECTS = 11000 
#             if bCustomRecords: 
#                 try:
#                     custom_records = pd.read_csv(f"{root_folder}/mimic3_annotations.csv",header = 0,index_col=0,dtype=str)
#                 except Exception as e:
#                     print(f'{e}, run "python src/lib/collate_datasets/create_mimic3_notes.py" in terminal to create the file')
#         elif MIMIC_NUM == "4":
#             MIMIC_PATH = 'mimic4wdb/0.1.0/'
#             MAX_SUBJECTS = 100

#         config = {
#             "paths": {
#                 "local": {
#                  "root_folder": Path(root_folder),
#                 },
#             },
#             "mimic_info": {
#                 "mimic_num": MIMIC_NUM,
#                 "mimic_path": MIMIC_PATH,
#             },
#             "version_num":VERSION_NUM,
#             "min_minutes" : MIN_MINUTES,
#             "num_subjects": MAX_SUBJECTS,
#             "categories_of_interest": [],
#             "custom_records": custom_records,
#             "ethnicity_extract": False,
#             'required_signals': {
#                     'ecg_labels' : ecg_labels,
#                     'ppg_labels' : ppg_labels,
#                     'abp_labels' : abp_labels,
#                     'custom_label': -1,
#             },
#         }

#         collator = MimicCollator(config,verbose=True)
#         collator.collate_dataset(load_waveforms=True)