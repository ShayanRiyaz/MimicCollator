# Mimic Collator
**MIMIC** (**M**edical **I**nformation **M**art for **I**ntensive **C**are) are large, freely-available database comprising deidentified health-related data associated with over patients who stayed in critical care units of the Beth Israel Deaconess Medical Center. This repository contains a parallelized, configurable Python class for extracting, filtering, and packaging MIMIC-III/IV waveform data (ECG, PPG, ABP) into HDF5. 

## Motivation & advantages for general use

- **Modularity**  
  - Single API supports both MIMIC-III and MIMIC-IV.  
  - Clear separation of concerns: path computation, cache lookup, signal checking, waveform I/O, batching.
- **Configurability**  
  - YAML-driven parameters: which signals, durations, subject counts, version numbers, etc.  
- **Scalability**  
  - Multi-threaded processing (via `ThreadPoolExecutor`) and HDF5 batching for large cohorts.
- **Reusability**  
  - Helper methods (`compute_paths`, `signal_requirements_check`, `build_subject_data`) can be adapted for other waveform repositories.
- **Ease of integration**  
  - Simply import `MimicCollator`, point it at your MIMIC data directory, and configure via YAML.

## Quick Start

1. **Clone the repository**  
```bash
git clone https://github.com/ShayanRiyaz/MimicCollator
cd MIMIC-Collator
```

2. **Set up environment and activate**
```bash
conda env create -f environment.yml
conda activate mimic_collator
```

3) **Modify ```config_mimic.yml``` if needed**
	-	Open config_mimic.yml and adjust parameters (e.g. required_signals, START_MINUTES, END_MINUTES, num_subjects).
```YAML
version_num: 1
min_minutes: 20
ethnicity_extract: False

subject_properties:
  custom_label: 0

required_signals:
  ecg_labels: ['ii']
  ppg_labels: ['pleth']
  abp_labels: ['abp', 'art']

mimic_variants:
  "3":
    mimic_path: 'mimic3wdb-matched/1.0/'
    max_subjects_available: 10999 #10999
    required_subjects: 10999 # Max Available 10999
    annotation_csv: 'mimic3_annotations.csv'
  "4":
    mimic_path: 'mimic4wdb/0.1.0/'
    max_subjects_available: 189
    required_subjects: 189 # 198
    annotation_csv: null
```

4) **Run**
```bash
python main.py
```
## Output .h5 Structure

```json
{
  "fix": {
    "subj_id": "<string>",
    "rec_id": "<string>",
    "files": "<string>",
    "af_status": "<int>",
    "subject_notes": "<any>"
  },
  "ppg": {
    "v":       "<number[]>",
    "fs":      "<int>",
    "method":  "<string>",
    "label":   "<string>"
  },
  "ekg": {
    "v":       "<number[]>",
    "fs":      "<int>",
    "method":  "<string>",
    "label":   "<string>"
  },
  "bp": {
    "v":       "<number[]>",
    "fs":      "<int>",
    "method":  "<string>",
    "label":   "<string>"
  }
}
```

## Key Concepts & Code Structure
-	collate_dataset() 
   Master workflow:
   ```python
   setup_paths → prepare record list → cache lookup → filter by signals → scan_waveforms → write batches
   ```
-	**Caching** (find_matching_records_info_file)
  -	Finds/validates JSON metadata caches.
  -	Picks the smallest cache ≥ desired_n, or largest < desired_n.
  - parallel record processing (extract_matching_record_info)
-	**Signal-presence tracking** (update_signal_presence)
  -	Maintains a one-hot history of which signals appear in each subject.
-	**Path resolution** (compute_paths)
  -	**Handles** both **local** filesystem and **PhysioNet** URL fall-back logic, and both MIMIC-III/IV layouts.
-	**Subject data builder** (build_subject_data)
  -	Packs metadata (fix) and waveforms (ppg, ekg, bp) into a nested dict.
-	**Parallel record processing** (scan_waveform_directory + _write_batch)
  -	Maps process_single_record across threads, batches results into HDF5.


## Configuration Options & Tuning
-	Window length: START_MINUTES, END_MINUTES
-	Signals: add/remove in required_signals
-	Batch size: MAX_SET_SIZE
-	Number of subjects: num_subjects, max_subjects_available
-	Verbosity: bVerbose

## Known Shortcomings & Future Work
-	Hard-coded signal support only ECG (lead II), PPG, ABP.
-	Rigid output schema: fixed keys (fix, ppg, ekg, bp)—no dynamic signal grouping.
-	Basic error handling: skips on exception, lacks retry or detailed logging.
-	~~No built-in filtering or QC: assumes cleaner signals or external preprocessing.~~ (Not the purpopse of this project)


## Contribution & License
- Special thanks to [Peter Charlton](https://github.com/peterhcharlton) whose work on Mimic III data collation  was of huge help for this project.
-	Issues: use GitHub Issues for bugs or feature requests.
-	License: MIT © Shayan Riyaz

## Citations
- Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.
- Gow, B., Pollard, T., Nathanson, L. A., Johnson, A., Moody, B., Fernandes, C., Greenbaum, N., Waks, J. W., Eslami, P., Carbonati, T., Chaudhari, A., Herbst, E., Moukheiber, D., Berkowitz, S., Mark, R., & Horng, S. (2023). MIMIC-IV-ECG: Diagnostic Electrocardiogram Matched Subset (version 1.0). PhysioNet. https://doi.org/10.13026/4nqg-sb35.
- Johnson, A. E. W., Pollard, T. J., Shen, L., Lehman, L. H., Feng, M., Ghassemi, M., Moody, B., Szolovits, P., Celi, L. A., & Mark, R. G. (2016). MIMIC-III, a freely accessible critical care database. Scientific Data, 3, 160035.
- Moody, B., Hao, S., Gow, B., Pollard, T., Zong, W., & Mark, R. (2022). MIMIC-IV Waveform Database (version 0.1.0). PhysioNet. https://doi.org/10.13026/a2mw-f949.
- Moody, B., Moody, G., Villarroel, M., Clifford, G. D., & Silva, I. (2020). MIMIC-III Waveform Database Matched Subset (version 1.0). PhysioNet. https://doi.org/10.13026/c2294b.
- S. K. Bashar, E. Ding, A. J. Walkey, D. D. McManus and K. H. Chon, "Noise Detection in Electrocardiogram Signals for Intensive Care Unit Patients," in IEEE Access, vol. 7, pp. 88357-88368, 2019, doi: 10.1109/ACCESS.2019.2926199.


