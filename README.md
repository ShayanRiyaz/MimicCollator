# MimicCollator
An independent MimicCollator class repository


## Steps to run:

### 1) Set up environment
```console
conda env create -f environment.yml
conda activate mimic_collator
```

### 2) Modify ```config_mimic.yml``` if needed
```YAML
version_num: 1
min_minutes: 30
ethnicity_extract: false

required_signals:
  ecg_labels: ['II']
  ppg_labels: ['Pleth']
  abp_labels: ['abp', 'art']
  custom_label: 0

mimic_variants:
  "3":
    mimic_path: 'mimic3wdb-matched/1.0/'
    max_subjects: 11000 # MAX 10999 available 
    annotation_csv: 'mimic3_annotations.csv'
  "4":
    mimic_path: 'mimic4wdb/0.1.0/'
    max_subjects: 100 # MAX 198 available
    annotation_csv: null
```

### 3) Run
```console
python main.py
```
