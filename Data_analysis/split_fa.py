with open('input.fa', 'r') as f:
    records = f.read().split('>')[1:] # split the file into records
    for i, record in enumerate(records):
        # create a new file for each record
        with open(f'{i+1}.fa', 'w') as output:
            output.write(f'>{record}') # write the record to the new file

