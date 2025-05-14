import pandas as pd
import os


mimic3_custom_dataset = {
    'label_map' : {
    0: 'normal sinus rhythm',
    1: 'atrial fibrillation',
}

}
mimic3_custom_dataset['filenames'] = ['p000608-2167-03-09-11-54','p000776-2184-04-30-15-16', 'p000946-2120-05-14-08-08', 'p004490-2151-01-07-12-36', 'p004829-2103-08-30-21-52', 'p075796-2198-07-25-23-40', 'p009526-2113-11-17-02-12', 'p010391-2183-12-25-10-15', 'p013072-2194-01-22-16-13', 'p013136-2133-11-09-16-58', 'p014079-2182-09-24-13-41', 'p015852-2148-05-03-18-39', 'p016684-2188-01-29-00-06', 'p017344-2169-07-17-17-32', 'p019608-2125-02-05-04-57', 'p022954-2136-02-29-17-52', 'p023824-2182-11-27-14-22', 'p025117-2202-03-15-20-28', 'p026377-2111-11-17-16-46', 'p026964-2147-01-11-18-03', 'p029512-2188-02-27-18-10', 'p043613-2185-01-18-23-52', 'p050089-2157-08-23-16-37', 'p050384-2195-01-30-02-21', 'p055204-2132-06-30-09-34', 'p058932-2120-10-13-23-15', 'p062160-2153-10-03-14-49', 'p063039-2157-03-29-13-35', 'p063628-2176-07-02-20-38', 'p068956-2107-04-21-16-05', 'p069339-2133-12-09-21-14', 'p075371-2119-08-22-00-53', 'p077729-2120-08-31-01-03', 'p087275-2108-08-29-12-53', 'p079998-2101-10-21-21-31', 'p081349-2120-02-11-06-35', 'p085866-2178-03-20-17-11', 'p087675-2104-12-05-03-53', 'p089565-2174-05-12-00-07', 'p089964-2154-05-21-14-53', 'p092289-2183-03-17-23-12', 'p092846-2129-12-21-13-12', 'p094847-2112-02-12-19-56', 'p097547-2125-10-21-23-43', 'p099674-2105-06-13-00-07']
mimic3_custom_dataset['labels'] = [0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1]



def make_label_map(filenames,df_labels,label_map,bMakeFile=0):
    labels = [label_map.get(label, 'unknown') for label in df_labels]
    df = pd.DataFrame({'filenames': filenames,'label': labels})
    df['subject_id'] = (df['filenames'].str.split('-', n=1).str[0].str.lstrip('p'))
    df['rec_id']  = df['filenames'].str.partition('-')[2]

    if bMakeFile:
        root_folder = "data/raw/mimic3_data/"
        if not os.path.exists(root_folder):
            os.makedirs(root_folder)
        df.to_csv(os.path.join(root_folder,'mimic3_annotations.csv'))
    return df

if __name__ == "__main__":
    bMakeFile = 1
    df = make_label_map(mimic3_custom_dataset['filenames'],mimic3_custom_dataset['labels'],mimic3_custom_dataset['label_map'],bMakeFile=bMakeFile)

