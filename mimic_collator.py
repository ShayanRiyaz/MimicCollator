'''
MimicCollator.py

Utilities for collating MIMIC III and IV waveform records into matfile sets.

Usage:
    from mimic_collator import MimicCollator
    # Set up config from config_mimic.yml
    collator = MimicCollator(config, verbose=True,search_cache=True)
    collator.collate_dataset(load_waveforms=True)

Author: Shayan Riyaz
License: MIT

'''

from rich import print
import os,time,json, requests, wfdb, h5py
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterator, Optional, Tuple, Dict, Any
from datetime import datetime, timezone
from urllib.request import urlretrieve

session = requests.Session()

def write_subject_record(h5_file : h5py.File, subj_id: str, data: Dict) -> None:
    """
    Write a single subject’s data into the HDF5 file under a dedicated group.

    This creates (or reuses) a group at /subjects/{subj_id}, then for each
    top-level section in `data` (e.g. "fix", "ppg", "ekg", "bp"), it creates a
    subgroup and writes each key→value pair as a dataset. Strings become
    fixed-length UTF-8 byte arrays; numeric scalars use the platform-independent
    default dtype; and NumPy arrays are GZIP-compressed and chunked.

    Args:
        h5_file (h5py.File):
            An open HDF5 file handle (in write mode).
        subj_id (str):
            The identifier for this subject. Determines the group path
            "subjects/{subj_id}".
        data (Dict[str, Dict[str, Union[str, int, float, np.ndarray]]]):
            Nested dictionary mapping section names to key→value dicts. 
            - `data[section][k]` may be:
                • str (will be stored as UTF-8 bytes),
                • int/float (stored as a scalar dataset),
                • np.ndarray (stored with gzip compression + chunking).

    Returns:
        None

    Raises:
        TypeError:
            If any `data[section][k]` is not one of the supported types
            (str, int/float, or np.ndarray).

    Notes:
        • Creating many small datasets can be slower; consider grouping large
          arrays in fewer fields if performance becomes an issue.
        • All datasets are created anew—if the key already exists, h5py will
          error. You can use `require_group` but not `require_dataset`.
        • Compression and chunking settings are tuned for waveform arrays.
    """
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

def get_typecast(cfg, folder_location=list, cast_fn=str, default=None, verbose=False):
    """
    Retrieve a nested config value and cast it, with a default fallback.

    Args:
        cfg (dict): Configuration dictionary to traverse.
        folder_location (List[str]): Sequence of keys indicating the nested path.
        cast_fn (callable): Function to apply to the raw value.
        default: Value to return if lookup or casting fails.
        verbose (bool): If True, print a warning on cast failure.

    Returns:
        The cast value at the specified path, or `default` if missing or invalid.
    """
    raw = cfg
    for k in folder_location:
        if isinstance(raw, dict):
            raw = raw.get(k, None)
        else:
            raw = None
        if raw is None:
            break
    try:
        return cast_fn(raw)
    except Exception:
        if verbose: print(f"Failed to cast {'/'.join(folder_location)}={raw}, defaulting to {default}")
        return default

class MimicCollator():
    """
    At the initial call for MIMIC III, this 
    """
    def __init__(self, config: dict,verbose: bool = False,search_cache: bool = True):
        """
        Initialize a MimicCollator for assembling waveform datasets.

        Args:
            config (dict)       : Nested configuration for paths, signal labels,
                sampling parameters, etc.
            verbose (bool)      : If True, print progress and debug info.
            search_cache (bool) : If True, attempt to load previously computed
                record‐lists from disk before re‐scanning.

        Attributes:
            Data Attributes:
            num_subjects (int)  : How many subjects to include.
            Constants:
                START_MINUTES (int), END_MINUTES (int): Window in minutes.
                START_SAMPLES (int), END_SAMPLES (int): Window in samples.
                MAX_SET_SIZE (int)                    : Number of subjets per set.
                MIN_FREQ    (float)                   : Minimum acceptable waveform sampling frequency

            Required channels:
                ecg_labels, ppg_labels, abp_labels (List[str])  : required ecg,ppg, abp signals
            Paths:
                root_dir(Path)                  : copied from config
                downloaded_files_dir (Path)     : copied from config 
                out_dir (Path)     : Where to save HDF5 files.
        """
        self.bVerbose = verbose
        self.search_cache = search_cache

        # Meta Data
        self.version_num = get_typecast(config, ["version_num"],cast_fn=int, default=0)
        self.ethnicity_extract = get_typecast(config,["ethnicity_extract"],cast_fn=str, default=False,verbose=self.bVerbose)   
        self.mimic_num = get_typecast(config, ["mimic_info","mimic_num"], cast_fn=str, default="3",verbose=self.bVerbose)        
        self.custom_label = get_typecast(config,["subject_properties","custom_label"], cast_fn=int, default=-1,verbose=self.bVerbose)    
        self.custom_records = get_typecast(config, ["custom_records"], cast_fn=lambda v: pd.DataFrame(v), default=pd.DataFrame())
        
        # Subject/sample parameters
        self.num_subjects = self.remaining_subjects = get_typecast(config, ["num_subjects"], cast_fn=int, default=50)
        self.max_subjects_available  = get_typecast(config, ["mimic_info","mimic_num","max_subjects_available"], cast_fn=int, default=100)

        # File paths
        self.root_dir  = get_typecast(config, ["paths","local","root_folder"],cast_fn=Path, default=os.getcwd())
        self.mimic_path  = get_typecast(config, ["mimic_info","mimic_path"],cast_fn=str, default="mimic4wdb/0.1.0/")
        self.downloaded_files_dir = self.root_dir / "downloaded_files"
        self.out_dir = self.root_dir / f'dataset_v{self.version_num}'
        if self.mimic_num == "4":
            self.mimic_matched_path = 'https://physionet.org/files/mimic-iv-ecg/1.0/'
        
        ## Constants
        self.MAX_SET_SIZE = get_typecast(config,[ "max_set_size"], cast_fn=int, default=50)
        self.MIN_FREQ = get_typecast(config, ["min_freq"], cast_fn=float, default=60.0)
        self.START_MINUTES = 5
        self.END_MINUTES = get_typecast(config, ["min_minutes"], cast_fn=int, default=20)
        self.START_SAMPLES = int(60*self.START_MINUTES*self.MIN_FREQ)
        self.END_SAMPLES = int(( self.END_MINUTES*60*self.MIN_FREQ)+ self.START_SAMPLES)
        
        # Signal Labels
        self.ecg_labels = [sig.lower() for sig in get_typecast(config, ["required_signals","ecg_labels"],cast_fn=list, default=[])]
        self.ppg_labels = [sig.lower() for sig in get_typecast(config, ["required_signals","ppg_labels"],cast_fn=list, default=[])]
        self.abp_labels = [sig.lower() for sig in  get_typecast(config, ["required_signals","abp_labels"],cast_fn=list, default=[])]
        # if (self.mimic_num == "3") and self.abp_labels:
        #     self.abp_labels = ['abp']

        self.all_required_signals = (self.ecg_labels or []) + (self.ppg_labels or []) + (self.abp_labels or [])
        self.required_signal_string = 'reqsigs_' + '_'.join(self.all_required_signals) 
        self.all_required_signals = self.all_required_signals or None
        
        # check all signals if present
        self.signal_presence_dicts = []  # One dict per subject
        self.all_signals_seen = set()    # Expanding set of signal labels

        if self.bVerbose:
            print(f'Minimum time: {self.END_MINUTES} minutes'
                  f'Required Signals: {self.required_signal_string}')


    def collate_dataset(self, load_waveforms=True):
        """
        Main Collate Dataset function. Steps:
        1) Initialize paths
        2) Build the list of candidate records (MIMIC-III or IV)
        3) Check for a cached mapping of “best” records
            a) If found and valid, load it
            b) Otherwise, regenerate via extract_matching_record_info
        4) Filter down to only those with all required signals
        5) Optionally scan/download waveform files and save
        """
        if self.bVerbose: print(f"Downloading and collating MIMIC {self.mimic_num} matched waveform subset")
        self.setup_paths()

        # Prepare Available Records
        if self.mimic_num == "3":
            match_records = self.prepare_mimic3_record_list()
        elif self.mimic_num == "4":
            match_records = self.prepare_mimic4_record_list()

        # Update available record length
        if self.max_subjects_available > match_records.shape[0]:
            if self.bVerbose: print(f'Not enough subjects available: Setting num_subjects from {self.num_subjects} to {match_records.shape[0]}')
            self.max_subjects_available = match_records.shape[0]

         # Check for existing records mathching length and metadata_requirement
        if self.search_cache:
            csv_path, existing_n  = self.find_matching_records_info_file(self.num_subjects)
        else:
            csv_path = []
            existing_n = None

        # Load Records if existing cache found
        if (existing_n is not None) and (self.max_subjects_available <= existing_n) and csv_path:
            df = pd.read_csv(csv_path)
            if self.bVerbose: print(f"Loaded {existing_n} from {csv_path.name}")
        else:
            # Create new matching records list
            if self.bVerbose: print("No best records CSV found; extracting new mapping")
            df = self.extract_matching_record_info(match_records)

        # Filter according to all_required_signals and num_subjects
        if self.all_required_signals is not None:
            df = df[df[self.required_signal_string] == 1]  
            if df.shape[0] < self.num_subjects:
                self.num_subjects = df.shape[0]
            df = df.head(self.num_subjects).reset_index(drop=True)
            if self.bVerbose: print(f"Loaded {existing_n}→sampled to {self.num_subjects} from df")
                
        # Load records from selected dataframe and save as h5
        if load_waveforms:
            self.scan_waveform_directory(match_records,record_info = df)
        return 


    def setup_paths(self) -> None:
        # Create root dir
        if not self.root_dir.exists():
            if self.bVerbose: print(f"Root folder does not exist: {self.root_dir}, setting up directories")
            self.root_dir.mkdir(parents=True, exist_ok=True)
      
        # Create Dowloaded files dir
        if not self.downloaded_files_dir.exists():
            if self.bVerbose: print(f"downloaded_files folder does not exist: {self.downloaded_files_dir}, setting up directories")
            self.downloaded_files_dir.mkdir(parents=True, exist_ok=True)
        
        # Create output files dir
        if not self.out_dir.exists():
            if self.bVerbose: print(f"Output folder does not exist: {self.out_dir}, setting up directories")
            self.out_dir.mkdir(parents=True,exist_ok=True)
        return



    def prepare_mimic3_record_list(self) -> Optional[pd.DataFrame]:
        """
            loads MIMIC III matched waveform repository record list from physionet, and custom labels if given .

            Outputs:
                df (DataFrame) : DataFrame of mimic3 subject id #'s.
        """

        df = {}
        recs = pd.read_csv(f"https://physionet.org/files/{self.mimic_path}/RECORDS-waveforms",header=None, names=["path"])
        df = pd.DataFrame({
            'matching_records':recs["path"].str.rpartition('/')[0] + '/',
            'subject_id': recs["path"].str.split("/", expand=True)[1].str.lstrip("p"),
            'filename': recs["path"].str.split("/", expand=True)[2]})
        
        df = df.drop_duplicates(subset='subject_id', keep='first')
        if not self.custom_records.empty: 
            ids_in_df2 = set(self.custom_records["subject_id"])
            df['in_df2'] = df['subject_id'].isin(ids_in_df2)
            any_matches = df['in_df2'].any()
            df['notes'] = str()
            if any_matches:
                label_map = self.custom_records.set_index('subject_id')['label'].to_dict()
                df['notes'] = df['subject_id'].map(label_map).fillna("").astype(str)
        df.reset_index(inplace=True,drop=True)
        return df
    




    def prepare_mimic4_record_list(self) -> Optional[pd.DataFrame]:
        """
            loads MIMIC IV waveform, and MIMIC IV ECG repository record list from physionet.
            Adds machine_measrument.csv notes from ECG dataset subjects to mactching subject
            ids in waveform dataset. 

            Outputs:
                df

            Notes: Current there is no match between MIMIC IV waveform and ECG repository
            subjects, so there won't be any subject notes.
        """
        df = {}
        dtype_cols = {i: str for i in range(16,22)}
        machine_measurments_path = self.downloaded_files_dir / 'machine_measurements.csv'
        url      = f"{self.mimic_matched_path}/machine_measurements.csv"

        if not machine_measurments_path.exists():
            urlretrieve(url, machine_measurments_path)   # C-optimized download
        machine_measurements_db = pd.read_csv(machine_measurments_path, dtype=dtype_cols)
        recs = pd.read_csv(f"https://physionet.org/files/{self.mimic_path}/RECORDS",header=None, names=["path"])

        df = pd.DataFrame({
            'matching_records':recs["path"],
            'subject_id':recs["path"].str.split("/", expand=True)[2].str.lstrip("p")
            })

        record_cols = [f"report_{i}" for i in range(18)]  
        def concat_records(group):
            data = group[record_cols].values.flatten()
            kept = [str(x) for x in data if pd.notna(x) and str(x).strip()!='']
            return ";".join(kept)
        
        ids_in_df2 = set(machine_measurements_db['subject_id'])
        df['in_df2'] = df['subject_id'].isin(ids_in_df2)
        any_matches = df['in_df2'].any()
        df["notes"] = ""
        if any_matches:
            mapping = (machine_measurements_db.groupby("subject_id").apply(concat_records).rename("all_records"))
            df = df[df["subject_id"].isin(mapping.index)].copy()
            df["all_records"] = df["subject_id"].map(mapping)
        return df
    



    def find_matching_records_info_file(self, desired_n: int) -> Tuple[Optional[Path], int]:
        """
        Find a metadata JSON cache meeting our requirements and return:

        Input:
            desired_n (int): Minimum number of subjects required.

        Output:
            Tuple[Optional[Path], int]: 
                - Path to the matching CSV file 
                - Number of subjects in that cache
                Returns (None, 0) if no valid cache is found.

        Notes:
            - Looks for `mimic{mimic_num}_*min_best_records_metadata.json` under root_dir.
            - Prefers the smallest cache ≥ desired_n; otherwise the largest < desired_n.
            - Verifies the corresponding CSV (replacing `_metadata.json` with `.csv`) exists.
        """
         
        pattern = Path(self.root_dir) / f"mimic{self.mimic_num}_*min_best_records_metadata.json"
        meta_files = sorted(pattern.parent.glob(pattern.name))

        supersets, subsets = [], []

        for meta_path in meta_files:
            try:
                with open(meta_path, 'r') as fp:
                    meta = json.load(fp)
                    
                subjects_available = meta.get("subjects_returned", 0)
                total_minutes = meta.get("end_minutes", 0)
                all_required_signals = meta.get("all_required_signals", 0)
                is_match = (subjects_available >= desired_n and
                            total_minutes >= self.END_MINUTES and
                            all_required_signals == self.all_required_signals)

                target_list = supersets if is_match else subsets
                target_list.append((subjects_available, meta_path))
            except Exception as e:
                if self.bVerbose: print(f"Skipping {meta_path} due to error: {e}")

        if supersets:
            selected = min(supersets, key=lambda x: x[0])
        elif subsets:
            selected = max(subsets, key=lambda x: x[0])
        else:
            selected = (None, 0)
            
        if selected[0] is None:
            return None, 0

        # Get corresponding CSV path from metadata file name
        _, meta_file = selected
        csv_file = meta_file.with_name(meta_file.name.replace("_metadata.json", ".csv"))

        if not csv_file.exists():
            if self.bVerbose: print(f"Found valid metadata at {meta_file} but missing CSV: {csv_file}")
            return None, 0

        return csv_file, selected[0]

    def update_signal_presence(self, actual_signals: list) -> Optional[Dict]:
        """
        Update signal presence records with the current subject’s signals.

        Input:
            actual_signals (list[str]):
                List of signal names present in the current subject.

        Output:
            dict[str, int]:
                A one-hot mapping of every signal seen so far (1 if present, 0 if absent) 
                for the current subject.

        Notes:
            - Any new signals are added to `self.all_signals_seen`.
            - Existing records in `self.signal_presence_dicts` are backfilled with 0 for new signals.
            - The returned dict is appended to `self.signal_presence_dicts`.
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
        '''
            This function:
            Read the .hea header, stripping out any '~' segments.
        '''
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

        if self.all_required_signals is not None:
            # signals_present = [sig in actual_signals for sig in self.all_required_signals]
            # all_signals_present = all(signals_present)
            # 1) Collect positions for each group
            ecg_positions = [i for i, s in enumerate(actual_signals) if s in self.ecg_labels]
            ppg_positions = [i for i, s in enumerate(actual_signals) if s in self.ppg_labels]
            abp_positions = [i for i, s in enumerate(actual_signals) if s in self.abp_labels]

            # 2) At least one from each?
            all_ok = bool(ecg_positions) and bool(ppg_positions) and bool(abp_positions)

            # 3) If it fails, return no-go
            if not all_ok:
                return hdr, 0, None
        return hdr,1,actual_signals




    def _extract_one_record(self, idx, curr_record, records_cache, df):
        """
        """
        rec_id = records_cache.get(curr_record)
        notes = df.loc[df["matching_records"] == curr_record,"notes"].squeeze() or ""
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
            
            info = {
                "idx":      idx,
                "pn_dir":   subject_and_rec_dir,
                "rec_id":   str(rec_id),
                "file_id":  segs[0][0],
                "seg_len":  segs[0][1],
                "max_freq": 0,
                "notes":    notes,
                f"{self.required_signal_string}": 0,
                "actual_signals": actual_signals,
            }
            return info

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
            info = {
                "idx":      idx,
                "pn_dir":   subject_and_rec_dir,
                "rec_id":   str(rec_id),
                "file_id":  seg_name,
                "seg_len":  hdr3.sig_len,
                "max_freq": hdr3.fs,
                f"{self.required_signal_string}": 1,
                "notes":    notes,
                "actual_signals": actual_signals,
            }
            return info
        return None

    def extract_matching_record_info(self, df, max_workers=8):
        '''
        
        '''
        matching_records = df['matching_records']
        if self.mimic_num == "3":
            records_cache = df.set_index('matching_records').to_dict()['filename']            
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
        args = [(idx, rec, records_cache, df) for idx, rec in enumerate(record_list)]
        infos = []
        # bDone = False
        num_subjects_observed = 0
        subject_count_with_all_requirements = 0
        try:
            with ThreadPoolExecutor(max_workers=max_workers) as exe:
                futures = {
                    exe.submit(self._extract_one_record, *arg): arg
                    for arg in args
                }
                for fut in as_completed(futures):
                    info = fut.result()
                    if info is not None:
                        # info['signal_labels'] = self.update_signal_presence(info['actual_signals'])
                        if info[self.required_signal_string] == 1:
                            subject_count_with_all_requirements +=1
                        
                            if self.bVerbose:
                                print(
                                f"Matching Subjects Found: {subject_count_with_all_requirements}, "
                                f"idx: {info['idx']}, "
                                f"pn_dir: {info['pn_dir']}, "
                                f"rec_id: {info['rec_id']}, "
                                f"file_id: {info['file_id']}, "
                                f"{self.required_signal_string}: {info[self.required_signal_string]} "
                                f"notes: {info['notes']}"
                                f"actual_signals: {info['actual_signals']}"
                                )

                        infos.append(info)
                        num_subjects_observed +=1
                        if subject_count_with_all_requirements >= self.num_subjects:
                            break
        except KeyboardInterrupt:
            if self.bVerbose: print("Interrupted. Cancelling futures...")
            for f in futures:
                f.cancel()
            raise

        
        df = pd.DataFrame(infos)
        # Sort Index
        df = df.sort_values(by=['idx'], ascending=True)
        # 1) Gather every unique signal across all subjects
        all_signals = sorted(df['actual_signals'].explode().dropna().unique())
    
        # 2) For each signal, add a column with 1 if present, 0 if not
        for sig in all_signals:
            df[sig] = df['actual_signals'].apply(lambda lst: int(sig in (lst or [])))
    
        # 3) Drop the raw list column
        df = df.drop(columns=['actual_signals'])
        
        # Create matching records csv file
        out_path = self.root_dir / f"mimic{self.mimic_num}_{self.END_MINUTES}min_best_records.csv"
        df.to_csv(out_path, index=False)

        if self.bVerbose: print(f"Wrote new {out_path.name}")
        meta = {
            "generated_at":  datetime.now(timezone.utc).isoformat() + "Z",
            "timestamp_unix": time.time(),
            "start_minutes": self.START_MINUTES,
            "end_minutes":   self.END_MINUTES,
            "ppg_labels":    self.ppg_labels,
            "ecg_labels":    self.ecg_labels,
            "abp_labels":    self.abp_labels,
            "all_required_signals": self.all_required_signals,
            "mimic_num":     self.mimic_num,
            "subjects_requested":  self.num_subjects,
            "subjects_returned":   num_subjects_observed,
            "num_subjects_with_all_requirements": subject_count_with_all_requirements,
        }

        meta_path = self.root_dir / f"mimic{self.mimic_num}_{self.END_MINUTES}min_best_records_metadata.json"
        with open(meta_path, "w") as mf:
            json.dump(meta, mf, indent=2)
        mf.close()
        if self.bVerbose:
            print(f"Saved metadata → {meta_path}")
        return df
    
    def get_records_iterator(self,df: Dict[str, Any],
        record_info: Optional[pd.DataFrame]) -> Iterator[Tuple[int,str]]:
        """
            Yield (idx, record_path) in the correct order, respecting record_info if given.
        """
        matching = df["matching_records"]
        if record_info is not None:
            matching = matching.iloc[record_info["idx"]]
        if hasattr(matching, "items"):
            yield from matching.items()
        else:
            yield from enumerate(matching)





    def compute_paths(self, curr_record: str, curr_idx: int, record_info: Optional[pd.DataFrame]) -> Tuple[str, str, str, str]:
        """
        Compute filesystem or URL paths for a given record and return:

        Input:
            curr_record (str): Identifier for the record subdirectory.
            curr_idx (int): Row index in `record_info` DataFrame.
            record_info (Optional[pd.DataFrame]): Lookup table with columns ['idx', 'rec_id', 'file_id'].

        Output:
            Tuple[str, str, str, str]:
                pn_dir   — Local or remote directory path containing the waveform files.
                hea_file — Header filename (e.g., "<file_id>.hea").
                rec_id   — Record identifier string.
                file_id  — File identifier string.

        Notes:
            - With `record_info`, looks up `rec_id` and `file_id` by index and constructs paths accordingly.
            - If `record_info` is None, falls back to downloading the RECORDS file from PhysioNet and parsing its first entry.
            - Accounts for differences in MIMIC-III vs. MIMIC-IV directory layouts.
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
        """
        Build structured data for a single subject.

        Input:
            record (object): WFDB-like record with attributes:
                - sig_name (list[str])
                - p_signal (np.ndarray)
                - fs (int)
            subj_id (str): Subject identifier.
            rec_id (str): Record identifier.
            file_id (str): File identifier.
            notes (Any): Arbitrary notes for this subject.

        Output:
            Dict[str, Any]: Nested dict with keys:
                - "fix": metadata fields (subj_id, rec_id, files, af_status, subject_notes)
                - "ppg", "ekg", "bp": each a dict with:
                    - v      (np.ndarray): signal values
                    - fs     (int): sampling rate
                    - method (str): description of extraction
                    - label  (str): signal label

        Notes:
            - The ABP label is selected from `self.abp_labels` matching `record.sig_name`.
            - Current structure is fixed and does not adapt to input variations.
            - Future work will allow custom dict layouts based on required signals and metadata.
        """
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
        Process one record tuple and return structured subject data or skip on failure.

        Input:
            args (tuple):
                - idx (int): Index of the record in the master list.
                - curr_record (str): Record directory identifier.
                - df (pd.DataFrame): Metadata DataFrame with 'subject_id' and 'notes'.
                - record_info (pd.DataFrame): Lookup with ['idx', 'rec_id', 'file_id'].

        Output:
            dict[str, Any] or None:
                - A subject_data dict (via build_subject_data) on success.
                - None if signal requirements fail or an exception is raised.

        Notes:
            - Resolves paths with compute_paths().
            - Validates signals with signal_requirements_check(); skips if pass_check == 0.
            - Loads waveform channels using wfdb.rdrecord with START/END minutes.
            - Extracts subj_id and notes, then delegates to build_subject_data().
            - Exceptions are caught; logs a skip message if bVerbose is True.
        """
        idx, curr_record, df, record_info = args
        record_info = record_info.sort_values(by=['idx'], ascending=True)
        try:
            pn_dir, hea_file, rec_id, file_id = self.compute_paths(
                curr_record, idx, record_info)

            file_header,pass_check,sigs = self.signal_requirements_check(filename=hea_file.split('.')[0],pn_dir=pn_dir,header_type='file')
            seg_name = file_header.record_name

            if pass_check == 0:
                return None
            ppg_i = [i for i, sig in enumerate(sigs) if sig.lower() in self.ppg_labels]
            ecg_i = [i for i, sig in enumerate(sigs) if sig.lower() in self.ecg_labels]
            abp_l = next(l for l in self.abp_labels if l in sigs)
            abp_i = sigs.index(abp_l)
            channels = [ppg_i,ecg_i,abp_i]
            channels = [x for sub in channels for x in (sub if isinstance(sub, list) else [sub])]
            record = wfdb.rdrecord(seg_name, pn_dir=pn_dir, sampfrom=int(self.START_MINUTES*60*file_header.fs), sampto=int(self.END_MINUTES*60*file_header.fs), channels = channels, return_res=32)

            subj_id     = curr_record.split("/")[-2]
            subj_notes  = df.loc[df['subject_id'] == subj_id[1:]]['notes'].values
            data        = self.build_subject_data( record, subj_id, rec_id, file_id, subj_notes)
            return data

        except Exception as e:
            if self.bVerbose: print(f"Skipping {curr_record!r}: {e}")
            return None

    def _write_batch(self, batch: list, set_num: int):
        """
        Write one batch of subject_data dicts into a single HDF5 file.

        Input:
            batch (list[dict]): Subject data dicts as returned by process_single_record.
            set_num (int): Batch index (0-based), used to number the output file.

        Output:
            None (creates an HDF5 file at the specified output directory).

        Notes:
            - Ensures output folder `out_dir/subject_{num_subjects}_v{version_num}` exists.
            - Filename pattern: 
                `mimic{mimic_num}_data_{num_subjects}_{set_num+1}.h5`
            - Uses `write_subject_record(h5f, subj_id, data)` to serialize each record.
            - If `bVerbose` is True, prints a confirmation message with batch size.
        """
        # Ensure the output folder exists
        set_folder = self.out_dir / f"subject_{self.num_subjects}_v{self.version_num}"
        set_folder.mkdir(parents=True, exist_ok=True)

        # Build filename
        fname = f"mimic{self.mimic_num}_data_{self.num_subjects}_{set_num + 1}.h5"
        out_path = set_folder / fname
        # Write all subjects in this batch
        with h5py.File(out_path, "w") as h5f:
            for data in batch:
                subj_id = data['fix']['subj_id']
                write_subject_record(h5f, subj_id, data)

        if self.bVerbose:
            print(f"Saved batch #{set_num + 1} → {fname} ({len(batch)} subjects)")

    def scan_waveform_directory(self, df: Dict[str,Any], record_info: Optional[pd.DataFrame], max_workers:int = 8) -> list:
        """
        Scan waveform records in parallel and write subject data in HDF5 batches.

        Input:
            df (pd.DataFrame): Metadata DataFrame of matching records.
            record_info (Optional[pd.DataFrame]): Lookup DataFrame with rec_id/file_id.
            max_workers (int): Number of worker threads.

        Output:
            None (writes HDF5 files to disk in batches).

        Notes:
            - Uses get_records_iterator() to generate (idx, record) tuples.
            - Runs process_single_record() in a ThreadPoolExecutor.
            - Accumulates up to MAX_SET_SIZE subjects per batch and calls _write_batch().
            - Stops when the desired number of subjects is processed.
            - Ensures the output directory exists.
            - Prints progress logs if bVerbose is enabled.
        """

        subject_batch      = []
        set_num           = 0
        total_subjects_added = 0
        num_subjects_added_in_batch = 0
        num_remaining_subjects = self.num_subjects
        arg_iter = [
            (idx, rec, df, record_info)
            for idx, rec in self.get_records_iterator(df, record_info)
        ]
        set_folder = self.out_dir / f'subject_{self.num_subjects}_v{self.version_num}'
        if not set_folder.exists():
            set_folder.mkdir(parents = True,exist_ok = True)

        with ThreadPoolExecutor(max_workers=max_workers) as exe:
            futures = {exe.submit(self.process_single_record, args): args for args in arg_iter}        
            for fut in as_completed(futures):
                data = fut.result()
                if data is None:
                    continue
                subject_batch.append(data)
                total_subjects_added +=1
                num_remaining_subjects -= 1
                
                num_subjects_added_in_batch += 1
                if (num_subjects_added_in_batch >= self.MAX_SET_SIZE) or (num_remaining_subjects == 0):
                    self._write_batch(subject_batch, set_num+1)
                    set_num += 1
                    num_subjects_added_in_batch = 0
                    subject_batch = []
                if self.bVerbose: print(f"Remainings Subject ({num_remaining_subjects}/{self.num_subjects-set_num*self.MAX_SET_SIZE}) | {num_subjects_added_in_batch} Subjects Added in Set #{set_num}")
                
                if total_subjects_added >= self.num_subjects:
                    break
        return