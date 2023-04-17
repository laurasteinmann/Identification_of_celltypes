from copy import copy
import numpy as np

def sort_matrix(df): #Sorts matrix: diagonal contains highest values based on column names
    sorting_steps = len(df.columns.values) #Sorting matrix for every existing column
    sorted_columns = []  # Create Filter for Columns
    sorted_rownames = []    #Create Filter for Rows
    to_sort = copy(df)
    for step in range(sorting_steps):
       if step < sorting_steps-1:#Sorting step for first sorting
           column_name = to_sort.columns.values[0]
           sorted = to_sort.sort_values(by = column_name, ascending = False)
           sorted_columns.append(sorted.columns.values[0])
           sorted_rownames.append(sorted.index.values[0])
           to_sort = sorted.iloc[1:,1:]
       else: #Last sorting step to append all of the rest rows step == sorting_steps-1
           column_name = to_sort.columns.values[0]
           sorted = to_sort.sort_values(by = column_name, ascending = False)
           sorted_columns.append(sorted.columns.values[0])
           sorted_rows = sorted.index.values
           for j in range(len(sorted.index.values)):
               sorted_rownames.append(sorted_rows[j])

    df = df[sorted_columns] #Sort dataframe by sorted columns
    df = df.reindex(sorted_rownames) #Sort dataframe by sorted rows
    return df

def normalizing_matrix(df):
    celltypes = len(df.columns)
    for celltype in range(celltypes):
        row = np.array(df.iloc[:,celltype])
        normalized = row/sum(row) *100
        df.iloc[:,celltype] = normalized
    return df