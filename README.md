# 🧬 WES Variant Interpretation Agent
### An AI-Powered Agent for Whole Exome Sequencing (WES) Variant Prioritization and Interpretation

**From Raw Variant Data to Evidence-Based Genomic Interpretation**

An intelligent, tool-using agent designed to automate and streamline the interpretation of Whole Exome Sequencing (WES) variants by integrating genomic annotation tools, population databases, pathogenicity prediction models, and Large Language Models (LLMs).

---

## 🚀 Overview

The **WES Variant Interpretation Agent** is a Python-based bioinformatics framework that processes Variant Call Format (VCF) files and generates structured, evidence-supported variant interpretation reports.

Instead of relying on a single prediction model, the system coordinates multiple bioinformatics tools and genomic knowledge bases to collect, evaluate, and prioritize variant-level evidence.

The architecture combines a rule-based Master Agent with an LLM-powered Clinical Impact Engine to produce interpretable and traceable results.

### Key Capabilities

- **Variant Quality Control (QC):** Filters low-quality variants based on QUAL, read depth (DP), FILTER, genotype (GT), and allele depth (AD).
- **Functional Annotation:** Retrieves gene-level and transcript-level consequences using Ensembl Variant Effect Predictor (VEP).
- **Clinical Evidence Retrieval:** Integrates ClinVar for reported clinical classifications.
- **Population Frequency Analysis:** Retrieves allele frequencies from gnomAD.
- **Pathogenicity Prediction:** Integrates precomputed AlphaMissense scores for missense variants.
- **Intelligent Tool Selection:** Uses conditional agent logic to determine which tools should be executed for each variant.
- **Variant Prioritization:** Calculates transparent heuristic priority scores based on collected evidence.
- **LLM-Assisted Interpretation:** Uses Gemini or OpenAI-compatible models to generate structured explanations of genomic evidence.
- **Automated Reporting:** Produces comprehensive JSON and human-readable Markdown reports.

---

## 🏗️ System Architecture

The framework follows a modular agent-based architecture:

```text
                WES VCF INPUT
                      |
                      v
              MASTER AGENT
                      |
          +-----------+-----------+
          |           |           |
          v           v           v
        QC Tool     VEP Tool   Evidence Tools
          |           |           |
          |           |      +----+---------+
          |           |      |      |       |
          |           |      v      v       v
          |           |   ClinVar gnomAD AlphaMissense
          |           |      |      |       |
          +-----------+------+------+
                      |
                      v
             EVIDENCE INTEGRATION
                      |
                      v
              VARIANT SCORING
                      |
                      v
            CLINICAL IMPACT ENGINE
                      |
               +------+------+
               |             |
               v             v
             Gemini       OpenAI
               |             |
               +------+------+
                      |
                      v
              FINAL REPORTS
                JSON + MD
```

### Agentic Workflow

The Master Agent operates through an iterative decision-making process:

1. **Perceive:** Reads the VCF file and extracts variant information.
2. **Plan:** Determines the appropriate tools for each variant.
3. **Act:** Executes annotation, evidence retrieval, and analysis.
4. **Observe:** Stores the results in the agent's evidence state.
5. **Reflect:** Evaluates missing evidence, conflicting results, and execution errors.
6. **Repeat or Finish:** Continues the analysis when necessary and generates the final reports.

Unlike a conventional fixed pipeline, the agent can adapt its execution strategy according to variant characteristics and previously collected evidence.

---

## 🔬 Integrated Bioinformatics Tools

| Tool / Database | Function |
|---|---|
| Ensembl VEP | Functional variant annotation |
| ClinVar | Clinical significance and review status |
| gnomAD v4 | Population allele frequency |
| AlphaMissense | Missense pathogenicity prediction |
| AlphaGenome | Optional non-coding variant analysis interface |
| pysam | VCF processing and genomic file handling |
| bcftools | Variant normalization (optional) |

### AI Models

| Model | Role |
|---|---|
| Gemini 2.5 Flash | LLM-based evidence interpretation |
| GPT-4o-mini | Alternative LLM interpretation engine |
| AlphaMissense | Deep learning-based missense scoring |
| Ensembl VEP | Rule-based functional consequence prediction |

---

## 🧠 Intelligent Variant Prioritization

The agent calculates a heuristic priority score using multiple evidence sources.

The scoring framework considers:

- VEP functional impact.
- ClinVar clinical classification.
- ClinVar review status and exact variant matching.
- gnomAD allele frequency.
- AlphaMissense pathogenicity scores.
- Conflicting or incomplete evidence.

Variants are ranked to support downstream expert review.

**Important:** The priority score is a project-specific heuristic and is NOT an ACMG/AMP classification or a validated clinical pathogenicity score.

---

## 📂 Project Structure

```text
wes_variant_agent/
│
├── requirements.txt
│
├── data/
│   ├── input/
│   │   └── sample.vcf
│   │
│   └── output/
│       ├── qc_report.json
│       ├── clinical_impact_report.json
│       └── clinical_impact_report.md
│
├── src/
│   ├── __init__.py
│   ├── config.py
│   ├── http_utils.py
│   ├── qc_tool.py
│   ├── vep_tool.py
│   ├── clinvar_tool.py
│   ├── gnomad_tool.py
│   ├── alphamissense_tool.py
│   ├── alphagenome_tool.py
│   ├── llm_client.py
│   ├── clinical_impact_engine.py
│   └── master_agent.py
│
└── tests/
    └── test_pipeline.py
```

---

## ⚙️ Installation

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/wes_variant_agent.git

cd wes_variant_agent
```

### 2. Create a virtual environment

```bash
python3 -m venv .venv

source .venv/bin/activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Configure the LLM provider

The framework supports three operating modes.

**Option A — Rule-based mode (No API required)**

```bash
export LLM_PROVIDER=none
```

**Option B — Google Gemini**

```bash
export LLM_PROVIDER=gemini
export GEMINI_API_KEY="YOUR_API_KEY"
```

**Option C — OpenAI**

```bash
export LLM_PROVIDER=openai
export OPENAI_API_KEY="YOUR_API_KEY"
```

Never commit API keys or sensitive genomic information to GitHub.

### 5. Optional: Configure AlphaMissense

Download the appropriate GRCh38 precomputed AlphaMissense dataset and its tabix index, subject to the dataset's licensing conditions.

```bash
export ALPHAMISSENSE_FILE="/path/to/AlphaMissense_hg38.tsv.gz"
```
## 🤖 Local AI Models — Ollama

The WES Variant Interpretation Agent supports a local LLM architecture through **[Ollama](https://ollama.com/)**, enabling privacy-conscious, on-device AI-powered genomic interpretation without requiring external LLM API keys.

Ollama provides a local inference environment for running open-weight language models directly on your workstation or server.

### Supported Model Architecture

| Component | Technology | Function |
|---|---|---|
| Local LLM Runtime | Ollama | Local model execution and inference |
| Primary LLM | Qwen / Llama / Gemma | Evidence-based variant interpretation |
| Functional Annotation | Ensembl VEP | Variant consequence annotation |
| Clinical Evidence | ClinVar | Clinical significance retrieval |
| Population Genetics | gnomAD | Allele frequency analysis |
| Deep Learning | AlphaMissense | Missense pathogenicity prediction |

### Recommended Local Models

Models can be selected according to available GPU memory, RAM, and computational resources.

- **Qwen:** General-purpose reasoning and structured output generation.
- **Llama:** Natural language interpretation and evidence summarization.
- **Gemma:** Lightweight local inference and structured explanations.

### Ollama Installation

Download and install Ollama:

https://ollama.com/download

Download a model:

```bash
ollama pull qwen2.5:7b
```

Start the model:

```bash
ollama run qwen2.5:7b
```

Ollama exposes a local REST API at:

```text
http://localhost:11434
```

### Local AI Configuration

The agent can communicate with Ollama through its local Chat API:

```python
import requests

OLLAMA_URL = "http://localhost:11434/api/chat"
MODEL = "qwen2.5:7b"

def local_llm(prompt):

    response = requests.post(
        OLLAMA_URL,
        json={
            "model": MODEL,
            "messages": [
                {
                    "role": "system",
                    "content": (
                        "You are a bioinformatics assistant "
                        "specialized in clinical genomics."
                    )
                },
                {
                    "role": "user",
                    "content": prompt
                }
            ],
            "stream": False,
            "format": "json"
        },
        timeout=180
    )

    response.raise_for_status()

    return response.json()["message"]["content"]
```

### Key Advantages

- Local AI inference without external LLM API keys.
- Greater control over genomic data privacy.
- Reduced dependency on third-party AI services.
- Flexible model selection.
- Integration with agent-based bioinformatics workflows.
- Potential deployment in private research and laboratory environments.

**Privacy Note:** Local inference keeps prompts on the local Ollama server, provided cloud models and external services are not used elsewhere in the workflow. Database annotation services may still require external network connections.
---

## ▶️ Usage

### Run the agent

```bash
python -m src.master_agent \
    --vcf data/input/sample.vcf \
    --top 10
```

### Run offline tests

```bash
python -m pytest -q
```

The project includes four offline tests covering QC, agent behavior, error handling, and LLM response parsing.

The offline tests use mock network tools and do not require API credentials.

---

## 📊 Output Reports

The system generates three main output files:

| Output | Description |
|---|---|
| `qc_report.json` | QC results, filtering decisions, and rejection reasons |
| `clinical_impact_report.json` | Ranked variants, collected evidence, flags, priority scores, and agent trace |
| `clinical_impact_report.md` | Human-readable interpretation report for prioritized variants |

Each prioritized variant can include:

- Chromosome and genomic position
- Reference and alternate alleles
- Gene symbol and transcript information
- HGVS nomenclature
- Functional consequences
- Population allele frequency
- ClinVar classification
- AlphaMissense score
- Priority score
- Supporting and conflicting evidence
- Confidence level
- Suggested further investigations

---

## 🛡️ Safety, Limitations & Responsible Use

This project is intended strictly for **educational and research purposes**.

- It is not a clinically validated diagnostic system.
- The priority score does not replace ACMG/AMP classification.
- LLM-generated interpretations are restricted to the evidence supplied to the model but still require expert verification.
- Missing database records do not automatically imply pathogenicity.
- Genomic reference assembly compatibility (GRCh38) must be verified.
- Patient identifiers and sensitive genomic information must not be transmitted to external services without appropriate authorization and privacy safeguards.
- All final variant interpretations must be reviewed by a qualified clinical geneticist.

---

## 🔭 Future Development

Potential extensions include:

- Integration of additional variant interpretation databases.
- ACMG/AMP evidence classification modules.
- Multi-agent architecture with specialized annotation and interpretation agents.
- LLM-based dynamic tool selection and planning.
- Integration of phenotype-driven prioritization using HPO terms.
- Automated FASTQ-to-VCF workflows.
- VCF-to-clinical-report automation.
- Advanced visualization of variant-level evidence.
- Integration with clinical genomics reporting workflows.

---

## 👨‍🔬 Author

**Mohammad Javad Naderi**

Bioinformatics Researcher | Biochemist | AI in Life Sciences

Research interests:
- AI-driven Bioinformatics
- Computational Genomics
- Next-Generation Sequencing (NGS)
- Clinical Variant Interpretation
- Deep Learning in Biology
- AI Agents for Healthcare

---

## 📜 License & Disclaimer

Please review the licensing requirements of individual databases and models, particularly AlphaMissense, before redistribution or commercial use.

This repository is an educational and research implementation. It must not be used as a standalone clinical diagnostic tool.

**Built for the intersection of Artificial Intelligence, Bioinformatics, and Clinical Genomics.**
