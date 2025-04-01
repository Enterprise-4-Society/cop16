import pandas as pd
from utils import load_data

###### LOAD DATA ####################################################################################################
# Load Jasanski ownership dataset
own_path = 'data/mining/open_database_mine_production/data/ownership.csv'

###### DATA CLEANING ################################################################################################
### Create the parent and subsite ids ###############################################################################
def clean_facility_ids(df):
    """Splits facility_id into parent and subsite facility ids."""
    df[['parent_facility_id', 'subsite_facility_id']] = df['facility_id'].str.split('.', expand=True)
    return df

### Clean ownership table ##########################################################################################
# 1) Modify strings for pattern identification

# Define replacement patterns
def apply_replacements(df):
    replacements = [
        ('%) and', '%),'),
        ('Medu Capital (15%) Minerales Y Productos Derivados S.A. (30%)', 'Medu Capital (15%), Minerales Y Productos Derivados S.A. (30%)'),
        ('Yancoal Australia Ltd (51%), Glencore Coal Pty Ltd. (49%).', 'Yancoal Australia Ltd (51%), Glencore Coal Pty Ltd. (49%)'),
        ('Sibanye-Stillwater (50%), Impala Platinum (50%(', 'Sibanye-Stillwater (50%), Impala Platinum (50%)'),
        ('Yancoal Australia Ltd (84.47%), Nippon Steel and Sumitomo Metal Australia Pty Limited (9.53%), Mitsubishi Materials (Australia) Pty Limited (6%).',
         'Yancoal Australia Ltd (84.47%), Nippon Steel and Sumitomo Metal Australia Pty Limited (9.53%), Mitsubishi Materials (Australia) Pty Limited (6%)'),
        ('Acacia Mining plc, 63.9% subsidiary of Barrick Gold Corporation', 'Acacia Mining plc (63.9%), subsidiary of Barrick Gold Corporation'),
        ('Bulyanhulu Gold Mine Limited (BGML) subsidiary of Acacia Mining plc (63.9%), subsidiary of Barrick Gold Corporation',
         'Bulyanhulu Gold Mine Limited (BGML) (63.9%), subsidiary of Acacia Mining plc subsidiary of Barrick Gold Corporation'),
        ('Pangea Minerals Limited (PML) subsidiary of Acacia Mining plc (63.9%), subsidiary of Barrick Gold Corporation',
         'Pangea Minerals Limited (PML) (63.9%), subsidiary of Acacia Mining plc subsidiary of Barrick Gold Corporation'),
        ('PT Nusa Halmahera Minerals (Newcrest Mining Limited 75%)', 'PT Nusa Halmahera Minerals (Newcrest Mining Limited) (75%)'),
        ('Freeport-McMoRan Inc. (72%), Sumitomo Metal Mining Arizona, Inc. (15%), SMM Morenci, Inc. (13%)',
         'Freeport-McMoRan Inc. (72%), Sumitomo Metal Mining Arizona Inc. (15%), SMM Morenci Inc. (13%)'),
        ('Lundin Mining Corp. (80%), Sumitomo Metal MiningCo., Ltd and Sumitomo Corporation (20%)',
         'Lundin Mining Corp. (80%), Sumitomo Metal MiningCo. Ltd and Sumitomo Corporation (20%)'),
        ('BHP (57.5%), Rio Tinto (30%), JECO Corporation consortium comprising Mitsubishi, JX Nippon Mining and Metals (10%), JECO2 Ltd (2.5%)',
         'BHP (57.5%), Rio Tinto (30%), JECO Corporation consortium comprising Mitsubishi JX Nippon Mining and Metals (10%), JECO2 Ltd (2.5%)'),
        ('Compañia de MinasBuenaventura S.A.A.', 'Compania de Minas Buenaventura S.A.A.'),
        ('Gecamines SA', 'Gecamines S.A.'),
        ("Newmont Mining Corporation", "Newmont Mining Corp."),
        ('After 2007: Goldcorp', 'Goldcorp')
    ]

    # Apply replacements using a loop
    for old, new in replacements:
        df['owners'] = df['owners'].str.replace(old, new, regex=False)
    return df

# 2) Keep only most recent owners
def keep_most_recent_owners(df):
    """Keeps only the most recent ownership records per facility."""
    # Step 1: Sort the DataFrame by 'facility_id' and 'year' to make sure years are ordered correctly
    df_sorted = df.sort_values(by=['facility_id', 'year'])

    # Step 2: Group by 'facility_id' and keep only the row with the most recent year in the main DataFrame
    # This will move all older years to the `previous_owners` DataFrame
    # `idxmax()` keeps the index of the maximum (most recent) year in each group
    most_recent_indices = df_sorted.groupby('facility_id')['year'].idxmax()

    # Create the main DataFrame with only the most recent owners
    return df_sorted.loc[most_recent_indices].reset_index(drop=True)

# 3) Allocate main-site owners when subsite owners is missing
def allocate_missing_owners(df, df_recent):
    """Allocates main-site owners when subsite owners are missing."""
    for _, row in df.iterrows():
        if row['facility_id'] not in df_recent['facility_id'].values:
            match = df_recent[(df_recent['parent_facility_id'] == row['parent_facility_id']) &
                              (df_recent['subsite_facility_id'] == '00')]
            if not match.empty:
                parent_row = match.iloc[0].copy()
                parent_row['facility_id'] = row['facility_id']
                parent_row['subsite_facility_id'] = row['subsite_facility_id']
                df_recent = pd.concat([df_recent, pd.DataFrame([parent_row])], ignore_index=True)
    return df_recent

# 4) Restructure the ownership table based on string patterns

def restructure_ownership(df):
    """Restructures the ownership table based on patterns in the 'owners' column."""
    pattern1 = r'^([\w\s\.,&-’\'()]+\(\d+\.?\d*%\))(, [\w\s\.,&-’\'()]+\(\d+\.?\d*%\))+$'
    pattern2 = r'^[\w\s\.,&-:’\'()]+\(\d+\.?\d*%\)$'
    pattern3 = r'^\w[\w\s\.,&-]+\(\d+\.?\d*%\), .+'

    # Pattern 1: Multiple owners with associated ownership in parenthesis and separated by a comma
    # Step 1: Identify rows where the format matches "owner1 (ownership1%), owner2 (ownership2%)"
    own_pattern1 = df[df['owners'].str.match(pattern1, na=False)]

    # Step 2: Split the 'owners' column for the identified rows
    # This split will separate based on commas followed by space, where each segment is "Owner (Ownership%)"
    split_owners = own_pattern1['owners'].str.split(r',\s*(?![^()]*\))', expand=True).stack().reset_index(level=1,
                                                                                                          drop=True)
    # Step 3: Create a new DataFrame to hold the split 'owners' column
    own_pattern1 = own_pattern1.drop(columns=['owners']).join(split_owners.rename('owners')).reset_index(drop=True)

    # Step 4: Extract ownership percentages using regex and assign them to a new 'ownership' column
    own_pattern1['ownership'] = own_pattern1['owners'].str.extract(r'(\d+\.?\d*)%')[0]

    # Step 5: Clean up the 'owners' column by removing ownership percentages, leaving only owner names
    own_pattern1['owners'] = own_pattern1['owners'].str.replace(r'\s*\(\d+\.?\d*%\)\s*', '', regex=True)

    # Pattern 2: Single owner with associated ownership in parenthesis
    # Step 1: Identify rows where the format matches "owner1 (ownership1%)"
    own_pattern2 = df[(~df['owners'].str.match(pattern1, na=False)) & (df['owners'].str.match(pattern2, na=False))]

    # Step 2: Extract ownership percentages using regex and assign them to a new 'ownership' column
    own_pattern2['ownership'] = own_pattern2['owners'].str.extract(r'\((\d+\.?\d*)%\)')[0]

    # Step 3: Remove ownership percentage values from the 'owners' column, leaving only the owner names
    own_pattern2['owners'] = own_pattern2['owners'].str.replace(r'\s*\(\d+\.?\d*%\)\s*', '', regex=True)

    # Pattern 3: Single owner with comment at the end of the string
    # Step 1: Identify rows where the format matches "owner1 (ownership1%), comment"
    own_pattern3 = df[(~df['owners'].str.match(pattern1, na=False)) & (~df['owners'].str.match(pattern2, na=False)) & (
        df['owners'].str.match(pattern3, na=False))]

    # Step 2: Extract ownership percentages using regex and assign them to a new 'ownership' column
    own_pattern3['ownership'] = own_pattern3['owners'].str.extract(r'\((\d+\.?\d*)%\)')[0]

    # Step 3: Extract comments from the owners column using regex
    # Comments start after a comma and space following the ownership percentage
    own_pattern3['comment'] = own_pattern3['owners'].str.extract(r', (.+)$')[0]

    # Step 4: Remove ownership percentage and comment values from the 'owners' column, leaving only the owner names
    own_pattern3['owners'] = own_pattern3['owners'].str.replace(r'\s*\(\d+\.?\d*%\), .+', '', regex=True).str.strip()

    return pd.concat([own_pattern1, own_pattern2, own_pattern3], ignore_index=True)


###### PROCESS DATA #################################################################################################

def process_mine_ownership(file_path=own_path):
    """Executes the full ownership data processing pipeline and returns the final DataFrame."""
    df = load_data(file_path)
    df = clean_facility_ids(df)
    df = apply_replacements(df)
    df_recent = keep_most_recent_owners(df)
    df_recent = allocate_missing_owners(df, df_recent)
    df_treated = restructure_ownership(df_recent)
    return df_treated

