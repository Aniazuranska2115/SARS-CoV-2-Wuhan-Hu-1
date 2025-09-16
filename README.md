# SARS-CoV-2 Variant Analysis Pipeline

Projekt zaliczeniowy do przedmiotu *Zarządzanie procesami analizy danych* (2024/2025).  
Celem projektu jest implementacja potoku analizy wariantów genetycznych wirusa SARS-CoV-2 z wykorzystaniem **Snakemake** i środowisk **Conda **.

---

##  Struktura katalogów
```
data/         # 
references/   # 
results/      # 
Snakefile     # 
config.yaml   # 
envs/         # 
```

---

## ⚙ Przygotowanie środowiska
Do uruchomienia potrzebne jest **Conda** (np. Miniconda/Miniforge) i **Snakemake**.

1. Zainstaluj Conda.
2. Stwórz bazowe środowisko tylko z `snakemake`:
   ```bash
   conda create -c conda-forge -c bioconda -n snakemake snakemake
   conda activate snakemake
   ```
3. Snakemake automatycznie utworzy dodatkowe środowiska Conda dla każdej reguły na podstawie plików `envs/*.yaml`.

---

## 📥Dane wejściowe
- Pobierz plik FASTQ odpowiadający próbce **SAMN33344176** i umieść w katalogu `data/`:
  ```
  data/SAMN33344176.fastq.gz
  ```
- Pobierz genom referencyjny **NC_045512.2** i umieść w katalogu `references/`:
  ```
  references/NC_045512.2.fa
  ```
- Przygotuj indeksy referencji:
  ```bash
  samtools faidx references/NC_045512.2.fa
  bwa index references/NC_045512.2.fa
  ```

---

## ▶ Uruchamianie workflow

### 1. Suchy bieg (sprawdzenie reguł, bez uruchamiania):
```bash
snakemake -n --use-conda
```

### 2. Uruchomienie analizy (4 rdzenie, per-rule conda):
```bash
snakemake --use-conda --conda-prefix ./.snakemake/conda -j 4
```

### 3. Generowanie raportu:
```bash
snakemake --report report.html
```

### 4. Generowanie grafu DAG (wymaga `graphviz`):
```bash
snakemake --dag | dot -Tpng > dag.png
```

---

##  Wyniki (w katalogu `results/`)
- `qc/` — kontrola jakości (FastQC, MultiQC)  
- `trim/` — przycięte sekwencje (fastp)  
- `*.sorted.bam` i indeksy BAM (BWA + Samtools)  
- `*.coverage.txt` — pokrycie  
- `*.vcf.gz` — wywołane warianty (bcftools)  
- `report.html` — raport końcowy Snakemake  
- `dag.png` — graf workflow  

---

## 📌 Uwagi
- Projekt korzysta z opcji **per-rule Conda envs**, dlatego każda reguła ma przypisany własny plik YAML w katalogu `envs/`.  
- Dzięki temu analiza jest w pełni **reprodukowalna** i nie wymaga instalacji wszystkich narzędzi w jednym środowisku.  
- Przy pierwszym uruchomieniu Snakemake pobierze wszystkie wymagane paczki (może to potrwać kilka minut).  
