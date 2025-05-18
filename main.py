import os
import yaml
import pandas as pd
from pathlib import Path
from mimic_collator import MimicCollator
from create_mimic3_notes import make_label_map, mimic3_custom_dataset

def load_yaml_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def build_config(base_cfg, mimic_num, root_folder, custom_records):
    return {
        "paths": {
            "local": {"root_folder": Path(root_folder)},
        },
        "mimic_info": {
            "mimic_num": mimic_num,
            "mimic_path": base_cfg["mimic_variants"][mimic_num]["mimic_path"],
            "max_subjects_available" : base_cfg["mimic_variants"][mimic_num]["max_subjects_available"]
        },
        "subject_properties": {
            "custom_label": base_cfg['subject_properties']["custom_label"]
        },
        "required_signals": {
            "ecg_labels" : base_cfg["required_signals"]["ecg_labels"],
            "ppg_labels" : base_cfg["required_signals"]["ppg_labels"],
            "abp_labels" : base_cfg["required_signals"]["abp_labels"],
        },
        "version_num": base_cfg["version_num"],
        "min_minutes": base_cfg["min_minutes"],
        "num_subjects": base_cfg["mimic_variants"][mimic_num]["required_subjects"],
        "categories_of_interest": [],
        "custom_records": custom_records,
        "ethnicity_extract": base_cfg["ethnicity_extract"],
        "required_signals": base_cfg["required_signals"]
    }

if __name__ == "__main__":
    base_cfg = load_yaml_config("config_mimic.yml")
 
    for mimic_num in base_cfg["mimic_variants"]:
        root_dir = Path(f'data/raw/mimic{mimic_num}_data/')
        

        custom_records = {}
        ann_csv = base_cfg["mimic_variants"][mimic_num]["annotation_csv"]

        if ann_csv:
            ann_path = Path(root_dir) / ann_csv
            make_label_map(mimic3_custom_dataset['filenames'],mimic3_custom_dataset['labels'],mimic3_custom_dataset['label_map'],bMakeFile=1)
            custom_records = pd.read_csv(ann_path, header=0, index_col=0, dtype=str)
        config = build_config(base_cfg, mimic_num, root_dir, custom_records)

        collator = MimicCollator(config, verbose=True,search_cache=True)
        collator.collate_dataset(load_waveforms=True)