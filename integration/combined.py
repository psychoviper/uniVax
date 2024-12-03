# Written by Shubham Prakash and Sufi during 01/10/24 to 27/11/24
# Purpose: This script will do four task.
# Purpose(cont): First will convert fasta file to csv file.
# Purpose(cont): Then using Algpred2 server, run allergicity test for certain number(INPUT BY USER) of sequence.
# Purpose(cont): Finally using Phobius server, will  detect presence of signal Peptides.
# Purpose(cont): Finally will filter out the sequence that succeed in all three above tests.

# Purpose(cont): THEN IT WILL RUN B CELL EPITOPES PREDICTION FOR ALL THE FILTERED SEQEUNCE
# Purpose(cont): Finally will run allergen, toxicity and antigenicity test.



#Input: fasta_file variable. I have named it input.fa. Feel free to change.
#Output: csv_file variable. I have named it sequence.csv. Feel free to change. It contains all sequence with their corresponding results.
#Output2: final_file variable. I have named it filtered_sequence.csv. Fell free to change. It contains only those sequence which  succeeded in allergen test, allergen test and signal P detection.
#Output3: 

# Libraries Needed
# lxml,bs4,pandas,

import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
from django.conf import settings
import os
######################################################################
# STEP-1: FASTA TO CSV
######################################################################
def parse_fasta(file):
    """Parse a FASTA file and return a dictionary with protein IDs and sequences."""
    # Initialize the dictionary to hold your data
    data = {
        'Protein ID': [],
        'Allergen Test':[],
        'Antigen Test':[],
        'Signal P':[],
        'Sequence': []
    }
    
    with open(file, 'r') as f:
        protein_id = ""
        sequence = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if protein_id:
                    data['Protein ID'].append(protein_id)
                    data['Sequence'].append(''.join(sequence))
                    data['Allergen Test'].append("Pending")
                    data['Antigen Test'].append("Pending")
                    data['Signal P'].append("Pending")
                protein_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
        if protein_id:
            data['Protein ID'].append(protein_id)
            data['Sequence'].append(''.join(sequence))
            data['Allergen Test'].append("Pending")
            data['Antigen Test'].append("Pending")
            data['Signal P'].append("Pending")
    return data
def fasta_to_csv(fasta_file, csv_file):
    """Convert FASTA to CSV using pandas."""
    data = parse_fasta(fasta_file)
    df = pd.DataFrame(data)
    df.to_csv(csv_file, index=False)

######################################################################
# STEP-2: ALLERGIC CHECK
######################################################################
def get_allergenicity_result(sequence,id,index):
    print("------------------------------------------------------")
    print(f"|   {index + 1}. Attempting: Allergen Test for: {id}  |")

    payload = {
        "name": "Job5",
        "seq": f">Protein\n{sequence}",  # Pass the sequence here
        "terminus": 4,
        "svm_th": 0.3,
    }
    url = "https://webs.iiitd.edu.in/raghava/algpred2/batch_action.php"
    
    r = requests.post(url, data=payload)
    
    if r.status_code == 200:
        soup = BeautifulSoup(r.text, 'html.parser')
        try:
            # Find the result in the table (you might need to adjust the index if it's different)
            # result = soup.find_all("td")[6].get_text()  # Assuming the 7th <td> contains the result
            result = soup.find_all("td")[6].get_text().strip()
        except IndexError:
            result = "ERROR:DATA EXTRACTION"
        print(f"|   {index + 1}. Success   : Allergen Test for: {id}  |")
    else:
        result = "ERROR:SERVER"
        print(f"|   {index + 1}. Failed   : Allergen Test for: {id}  |")
    return result
######################################################################
# STEP-3: Signal P Detection
######################################################################
def get_signalP_result(sequence, id, index):
    print("------------------------------------------------------")
    print(f"|      {index + 1}. Detecting: SignalP for: {id}      |")
    payload = {
        "protseq": f">{id}\n{sequence}",
        "format": 'nog',
    }
    url = "https://phobius.sbc.su.se/cgi-bin/predict.pl"
    
    r = requests.post(url, data=payload)
    
    if r.status_code == 200:
        soup = BeautifulSoup(r.text, 'html.parser')
        try:
            text = soup.find("pre").get_text()
            if "signal" in text.lower():
                result = "Found"
            else:
                result = "Not Found"
        except IndexError:
            result = "ERROR:DATA EXTRACTION"
        print(f"|      {index + 1}. Success  : SignalP for: {id}      |")
    else:
        result = "ERROR:SERVER"
        print(f"|       {index + 1}. Failed  : SignalP for: {id}      |")
    return result
######################################################################
# STEP-5: B Cell Epitope Detection and allergen, toxicity, antigenicity test on those epitopes.
######################################################################
def bcell(final_file):
    ######################################################################
    # FUNCTION TO SEND POST REQUEST TO ABCPred
    ######################################################################
    def get_bcell_epitof(sequence, id, index, d2, threshold):
        print("------------------------------------------------------")
        print(f"|   {index + 1}. Attempting: B Cell Epitope for: {id}  |")
        url = "https://webs.iiitd.edu.in/cgibin/abcpred/test1_main.pl"
        payload = {
            "SEQNAME": "",
            "SEQ": f"{sequence}",
            "Threshold": threshold,
            "window": 16,
            "filter": "on"
        }
        r = requests.post(url, data=payload)
        epitope_count = 0  # Initialize the count of epitopes
        if r.status_code == 200:
            soup = BeautifulSoup(r.text, 'lxml')
            try:
                i = 0
                for elements in soup.find_all('td', attrs={'width': '50%'}):
                    if i == 0:  # Skip the first irrelevant row
                        i += 1
                        continue
                    d2['Protein ID'].append(id)
                    d2['Epitope'].append(elements.get_text())
                    i += 1
                    epitope_count += 1
                score_elements = soup.find_all('td', attrs={'width': '20%'}) # For score calculation
                for i in range(3, len(score_elements), 2):  # Start at 3, step by 2
                    d2['Score'].append(score_elements[i].get_text(strip=True))
            except:
                d2['Protein ID'].append(id)
                d2['Epitope'].append("ERROR:Data Extraction")
                d2['Score'].append("NA")
            print(f"|   {index + 1}. Success   : B Cell Epitope for: {id}  |")
        else:
            print(f"|    {index + 1}. Failed  : B Cell Epitope for: {id}   |")
        return epitope_count  # Return the count of epitopes found

    ######################################################################
    # FUNCTION TO PRINT COUNT OF EPITOPS FOR SAKE OF RE-ITERATION IF NEEDED
    ######################################################################
    def print_summary(d3,threshold):
        print("------------------------------------------------------")
        print("\n\n------------------------------------------------------")
        print(f"|    Summary of Epitope Counts (Threshold: {threshold})    |")
        print("------------------------------------------------------")
        print(f"{'Protein ID':<26}{'Epitope Count':>26}")
        print("------------------------------------------------------")
        for i in range(len(d3['Protein ID'])):
            print(f"{d3['Protein ID'][i]:<26}{d3['Epitope Count'][i]:>26}")
        print("------------------------------------------------------")

    ######################################################################
    # FUNCTION TO RUN ALLERGIC TEST FOR A EPITOPE
    ######################################################################
    def get_allergenicity_result(sequence,id,index):
        print("------------------------------------------------------")
        print(f"|   {index + 1}. Attempting: Allergen Test for: {sequence}  |")

        payload = {
            "name": "Job5",
            "seq": f">Protein\n{sequence}",  # Pass the sequence here
            "terminus": 4,
            "svm_th": 0.3,
        }
        url = "https://webs.iiitd.edu.in/raghava/algpred2/batch_action.php"
        
        r = requests.post(url, data=payload)
        
        if r.status_code == 200:
            soup = BeautifulSoup(r.text, 'html.parser')
            try:
                # Find the result in the table (you might need to adjust the index if it's different)
                # result = soup.find_all("td")[6].get_text()  # Assuming the 7th <td> contains the result
                result = soup.find_all("td")[6].get_text().strip()
            except IndexError:
                result = "ERROR:DATA EXTRACTION"
            print(f"|   {index + 1}. Success   : Allergen Test for: {sequence}  |")
        else:
            result = "ERROR:SERVER"
            print(f"|   {index + 1}. Failed   : Allergen Test for: {sequence}  |")
        return result
    ######################################################################
    # FUNCTION TO RUN TOXICITY TEST FOR A EPITOPE
    ######################################################################
    def get_toxicity(output_file,num_sequences):
        ######################################################################
        # FUNCTION TO SEND POST REQUEST TO TOXINPRED SERVER
        ######################################################################
        def toxicity_result(sequence,id,index):
            url = "https://webs.iiitd.edu.in/raghava/toxinpred/multiple_test.php"
            cookies = {
                "PHPSESSID": "75s8sth03frii1hosaer2f39j1",  # Replace with the actual session ID if required
            }
            headers = {
                "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36",
                "Referer": "https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php",
            }
            data = {
                "seq": sequence,
                "method": "1",
                "eval": "10",
                "thval": "0.0",
                "field[]": ["4", "7", "9", "11", "13"],
            }
            print("------------------------------------------------------")
            print(f"|   {index + 1}. Attempting: Toxicity Test for: {sequence}  |")
            with requests.Session() as session:
                # Step 1: Submit the POST request
                response = session.post(url, data=data, headers=headers, cookies=cookies)
                
                # Check the response
                if response.status_code == 200:
                    # Step 2: Parse the meta-refresh URL
                    soup = BeautifulSoup(response.text, "html.parser")
                    meta_refresh = soup.find("meta", attrs={"http-equiv": "refresh"})
                    
                    if meta_refresh:
                        # Extract the URL from the content attribute
                        refresh_url = meta_refresh["content"].split("url=")[-1]
                        final_url = f"https://webs.iiitd.edu.in/raghava/toxinpred/{refresh_url}"
                        
                        # Step 3: Send a GET request to the redirected URL
                        final_response = session.get(final_url)
                        
                        if final_response.status_code == 200:
                            # print("Final Response Text:\n", final_response.text)
                            soup=BeautifulSoup(final_response.text,"html.parser")
                            result=soup.find_all("td", align="center")[3].get_text() #Our final result (index need to be checked routinely)
                            print(f"|   {index + 1}. Success   : Toxicity Test for: {sequence}  |")
                        else:
                            result = "ERROR:DATA EXTRACTION"
                            print(f"|   {index + 1}. Failed   : Toxicity Test for: {sequence}  |")

                    else:
                        result = "ERROR:NO REDIRECTION PAGE FOUND"
                        print(f"|   {index + 1}. Failed   : Toxicity Test for: {sequence}  |")

                else:
                    result = "ERROR: INITIAL REQUEST FAILED"
                    print(f"|   {index + 1}. Failed   : Toxicity Test for: {sequence}  |")
                return result
            
        ######################################################################
        # DRIVER CODE STARTS
        ######################################################################
        input_file=output_file

        df=pd.read_csv(input_file)
        for index, row in df.head(num_sequences).iterrows():
            sequence = row['Epitope']
            id=row['Protein ID']
            result = toxicity_result(sequence,id,index) 
            # Add the 'Toxicity Test' column with the fetched result
            df.at[index, 'Toxicity Test'] = result
        df.to_csv(input_file, index=False)
    ######################################################################
    # FUNCTION TO RUN ANTIGEN TEST FOR A EPITOPE (!!!PENDING!!!!)
    ######################################################################
    def get_antigen(output_file,num_sequences):
        ######################################################################
        # FUNCTION TO SEND POST REQUEST TO VAXIJEN SERVER
        ######################################################################
        def antigenicity_result(sequence,id, index):
            return "Pending"
        
        ######################################################################
        # DRIVER CODE STARTS
        ######################################################################
        input_file=output_file
        df=pd.read_csv(input_file)
        for index, row in df.head(num_sequences).iterrows():
            sequence = row['Epitope']
            id = row['Protein ID']
            result = antigenicity_result(sequence, id, index)
            df.at[index, 'Antigen Test'] = result
        df.to_csv(input_file, index=False)

    ######################################################################
    # MAIN CODE STARTS
    ######################################################################
    start_time = time.time()
    csv_file = "filtered_sequence.csv"  # INPUT CSV File. Update Accordingly.
    df = pd.read_csv(csv_file)
    d2 = {'Protein ID': [], 'Epitope': [], 'Score': []}  # New dictionary to store epitopes.
    d3 = {'Protein ID':[], 'Epitope Count':[]} # New Dictionary to store count of epitopes.

    # Ensure the number of sequences to process is within the available range
    num_sequences = int(input(f"|            Total sequences available: {len(df)}        |\n>How many sequences do you want to process for B Cell Epitope Detection? "))
    num_sequences = min(num_sequences, len(df))

    # Process only the specified number of sequences
    while True:
        threshold = float(input(f"> Enter threshold (default is 0.51): ") or 0.51)
        for index, row in df.head(num_sequences).iterrows():
            sequence = row['Sequence']
            id = row['Protein ID']
            epitope_count = get_bcell_epitof(sequence, id, index, d2, threshold)
            d3['Protein ID'].append(id)
            d3['Epitope Count'].append(epitope_count)
        print_summary(d3, threshold)
        action = input("> Do you want to repeat the process? (y/n): ").strip().lower() or 'n'

        if action != 'y':
            break
        else: #if wants to repeat the process then empty the dictionary containing epitopes and their count
            for key in d2:
                d2[key] = []
            for key in d3:
                d3[key] = []

            

    # Convert the updated d2 dictionary to df2 dataframe
    df2 = pd.DataFrame(data=d2) #contains epitopes
    output_file = "bCell.csv"  # OUTPUT CSV File. Update Accordingly.
    df2.to_csv(output_file, index=False)
    ######################################################################
    # NOW NEW STEPS FOR THREE TESTS OF THE EPITOPES THAT WE GOT!
    ######################################################################
    ######################################################################
    # ALLERGEN TEST
    ######################################################################
    df = pd.read_csv(output_file)
    num_sequences = int(input(f"|            Total Epitopes available: {len(df)}        |\n>How many Epitopes do you want to process for allergen and toxicity test? "))
    num_sequences = min(num_sequences, len(df))
    for index, row in df.head(num_sequences).iterrows():
        sequence = row['Epitope']  # Get the sequence
        id=row['Protein ID']
        result = get_allergenicity_result(sequence,id,index) 
        
        # Add the 'Allergen Test' column with the fetched result
        df.at[index, 'Allergen Test'] = result
    df.to_csv(output_file, index=False)
    ######################################################################
    # TOXICITY TEST
    ######################################################################
    get_toxicity(output_file,num_sequences)

    ######################################################################
    # ANTIGEN TEST
    ######################################################################
    get_antigen(output_file,num_sequences)

    ######################################################################
    # FILTERING EPITOPES WHICH SUCCEED IN ALLERGEN, ANTIGEN AND TOXICITY TEST
    ######################################################################
    # Read the original CSV file
    df = pd.read_csv(output_file)

    # Apply the filter conditions
    filtered_df = df[
        (df['Allergen Test'] == 'Non-Allergen') & 
        (df['Toxicity Test']=='Non-Toxin') &
        (df['Antigen Test'] == 'Pending')#change this after proper implementation of antigen-test (vaxijen)
    ]

    # Save the filtered data to a new CSV file
    output_csv = 'filtered_bCell.csv'  # OUTPUT CSV File. Change Accordingly.
    filtered_df.to_csv(output_csv, index=False)


    end_time = time.time()
    total_time = end_time - start_time
    print("------------------------------------------------------")
    print(f"| Success: Detecting B Cell Epitopes.                |")
    print(f"| Success: Allergen Test of Epitopes.                |")
    print(f"| Success: Toxicity Test of Epitopes.                |")
    print(f"| PENDING: Antigen Test of Epitopes.                 |")
    print(f"| Check {output_file} & {output_csv}                    |")       
    print(f"| Time Taken: {total_time:.2f} sec                              |")
    print("======================================================")

######################################################################
# MAINNNNNNNN CODE STARTS
######################################################################
csv_file = 'sequence.csv'  # OUTPUT CSV File. Change Accordingly.
fasta_file = 'input.fa'  # INPUT FASTA File. Change Accordingly.
wholeTimeStart=time.time()

def f_to_c(file_path):
    fasta_file = file_path  # INPUT FASTA File. Change Accordingly.
    csv_file = os.path.join(settings.MEDIA_ROOT,'sequence.csv') # OUTPUT CSV File. Change Accordingly.
    print("======================================================")
    print("|                 PROCESS INITIATED                  |")
    print("======================================================")
    print("|                                                    |")
    print("| Step1: FASTA to CSV Conversion                     |")
    print("| Step2: Allergenicity Test (AlgPred2)               |")
    print("| Step3: Signal Peptide Detection (Phobius)          |")
    print("| Step4: Filtering Sequence which succeed in.....    |")
    print("| .......in Allergen Test and Signal P Detection     |")
    print("|                                                    |")
    print("======================================================")
    print("|                  STEP-1 INTIATED                   |")
    print("======================================================")
    wholeTimeStart=time.time()
    start_time = time.time()
    fasta_to_csv(fasta_file, csv_file)
    end_time = time.time()
    total_time = end_time - start_time
    print(f"| Success: Converted {fasta_file}->{csv_file}  |")
    print(f"| Time Taken: {total_time:.2f} sec                               |")
    print("======================================================")
    print("|                  STEP-1 COMPLETED                  |")
    print("======================================================\n\n")

    print("======================================================")
    print("|                  STEP-2 INTIATED                   |")
    print("======================================================")
    df = pd.read_csv(csv_file)
    # num_sequences = int(input(f"|            Total sequences available: {len(df)}        |\n|     How many sequences do you want to process? "))
    # num_sequences = min(num_sequences, len(df))
    num_sequences = len(df)
    start_time = time.time()
    for index, row in df.head(num_sequences).iterrows():
        sequence = row['Sequence']
        id = row['Protein ID']
        result = get_allergenicity_result(sequence, id, index)
        df.at[index, 'Allergen Test'] = result
    df.to_csv(csv_file, index=False)
    end_time = time.time()
    total_time = end_time - start_time
    print("------------------------------------------------------")
    print(f"| Success: Done with Allergen Test.                  |")
    print(f"| Check {csv_file}                              |")       
    print(f"| Time Taken: {total_time:.2f} sec                              |")
    print("======================================================")
    print("|                  STEP-2 COMPLETED                  |")
    print("======================================================\n\n")

# if __name__ == '__main__':
    
    # print("======================================================")
    # print("|                  STEP-3 INTIATED                   |")
    # print("======================================================")
    # start_time = time.time()
    # # df = pd.read_csv(csv_file)
    # for index, row in df.head(num_sequences).iterrows():
    #     sequence = row['Sequence']
    #     id = row['Protein ID']
    #     result = get_signalP_result(sequence, id, index)
    #     df.at[index, 'Signal P'] = result
    # df.to_csv(csv_file, index=False)
    # end_time = time.time()
    # total_time = end_time - start_time
    # print("------------------------------------------------------")
    # print(f"| Success: Done with SignalP Detection.              |")
    # print(f"| Check {csv_file}                              |")       
    # print(f"| Time Taken: {total_time:.2f} sec                               |")
    # print("======================================================")
    # print("|                  STEP-3 COMPLETED                  |")
    # print("======================================================\n\n")
    # print("======================================================")
    # print("|                  STEP-4 INTIATED                   |")
    # print("======================================================")
    # start_time = time.time()
    # df = pd.read_csv(csv_file)
    # final_file='filtered_sequence.csv' #FINAL OUTPUTFILE AFTER STEP123
    # filtered_df = df[
    #     (df['Allergen Test'] == 'Non-Allergen') & 
    #     (df['Antigen Test'] == 'Pending') & #Change this after proper implementation of antigen-test (vaxijen)
    #     (df['Signal P'] == 'Found')
    # ]
    # filtered_df.to_csv(final_file, index=False)
    # end_time = time.time()
    # total_time = end_time - start_time
    # print("------------------------------------------------------")
    # print(f"| Success:Filtered sequence which are non-Allergen, have antigen and Signal P.|")
    # print(f"| Check {final_file}                              |")       
    # print(f"| Time Taken: {total_time:.2f} sec                               |")
    # print("======================================================")
    # print("|                  STEP-4 COMPLETED                  |")
    # print("======================================================\n\n")
    # print("======================================================")
    # print("|                  STEP-5 INTIATED                   |")
    # print("======================================================")
    # bcell(final_file)
    # print("======================================================")
    # print("|                  STEP-5 COMPLETED                  |")
    # print("======================================================\n\n")
    # wholeTimeEnd=time.time()
    # wholeTotalTime=wholeTimeEnd-wholeTimeStart
    # print(f"| All five steps successfully completed.            |")
    # print(f"| Check {csv_file} and {final_file}. Thank You!  |")       
    # print(f"| Total Time Taken: {wholeTotalTime:.2f} sec                        |")
    # print("======================================================")
    # print("|                PROCESS TERMINATED                  |")
    # print("======================================================\n\n")