# uniVax: Vaccine Design Software

**uniVax** is an integrated software platform designed to streamline the entire vaccine development process, covering key stages such as allergenicity prediction, antigenicity prediction, toxicity assessment, epitope prediction, and signal peptide detection. This tool combines over 50 different software tools into one cohesive workflow, reducing the need for users to jump between different platforms and tools. 

By leveraging **web scraping technologies** and **concurrency techniques**, **uniVax** significantly optimizes the vaccine design process, reducing what traditionally takes months into just hours or days.

The platform provides comprehensive control to the user while abstracting the complexity of the underlying processes, allowing for a seamless experience.
![Screenshot (138)](https://github.com/user-attachments/assets/c3563da0-fc9f-4497-8703-ae5478704d62)

### Key Features

1. **Allergenicity Prediction**: 
   - Predicts whether a given protein sequence is allergenic, based on well-established databases and prediction models.

2. **Antigenicity Prediction**: 
   - Utilizes state-of-the-art tools to predict the antigenicity of a protein sequence, ensuring that only the most relevant candidates are considered for vaccine development.

3. **Toxicity Assessment**:
   - Assesses the potential toxicity of a given sequence, helping to filter out unsafe candidates early in the process.

4. **Signal Peptide Detection**:
   - Identifies signal peptides within protein sequences that may play a role in directing the protein to the correct cellular location.

5. **Epitope Prediction**:
   - Predicts **B-cell epitopes**, **Helper T-cell epitopes**, **A1 Supertype** epitopes, and **B58 Supertype** epitopes to help target immune responses effectively.

6. **Filtration for Toxicity, Allergenicity, and Antigenicity**:
   - Filters epitope sequences based on toxicity, allergenicity, and antigenicity to ensure only the most promising candidates remain for vaccine design.

7. **Web Scraping & Concurrency**:
   - Uses advanced web scraping to gather necessary data and prediction results from multiple sources.
   - Implements concurrency to speed up processes, ensuring a highly optimized and fast workflow that saves valuable research time.
![Screenshot (139)](https://github.com/user-attachments/assets/575eac22-b7c4-434c-ac07-8ea46b0c3f70)

### Time Optimization
Previously, users had to visit multiple websites, input the same protein sequences repeatedly, and wait for long periods for results. **uniVax** automates and integrates these steps, cutting the process from months to mere hours or days. The software handles everything from allergenicity and antigenicity predictions to epitope linkage, all while saving the user from manually handling huge datasets.
![manual vs automation-1](https://github.com/user-attachments/assets/dcd1c2ed-1ee2-425d-add4-d3442e10f904)

### User Control and Workflow Abstraction

Despite the complexity of the underlying steps, **uniVax** is designed with a user-friendly interface. Necessary prompts and controls are provided, so users can easily track progress and make informed decisions throughout the process. This abstraction simplifies vaccine design for researchers, even those with limited technical backgrounds.
![Screenshot (137)](https://github.com/user-attachments/assets/a4622209-635b-45ca-85c8-2c14f896351c)

### Deployment
**uniVax** has been deployed and is available for use at the following URL:
[https://automation.sbhmprksh.in](https://automation.sbhmprksh.in)

### Installation and Setup

1. **Requirements**:
   - Python 3.8+
   - Required libraries: `Flask`, `requests`, `BeautifulSoup`, `concurrent.futures`, `NumPy`, `pandas`, and others (refer to `requirements.txt` for the complete list).
   
2. **Clone the Repository**:
   - Clone this repository to your local machine using:
     ```bash
     git clone https://github.com/your-repository/uniVax.git
     ```
   
3. **Install Dependencies**:
   - Navigate to the project directory and install the required dependencies:
     ```bash
     pip install -r requirements.txt
     ```

4. **Run the Application**:
   - Start the server with:
     ```bash
     python app.py
     ```

5. **Access the Application**:
   - Open your web browser and navigate to `http://localhost:5000` to access the **uniVax** interface.

### Contributing

We welcome contributions to **uniVax**! If you'd like to help improve the software, please fork the repository, create a branch for your feature or bugfix, and submit a pull request.

---

**Note**: Ensure your machine has internet access as the software utilizes web scraping for fetching necessary data.
