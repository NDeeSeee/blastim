# Глоссарий VEP: Consequences и поля VCF

## VCF — фиксированные поля

| Колонка | Описание | Пример |
|---------|----------|--------|
| CHROM | Хромосома | `chr17` |
| POS | Позиция (1-based, GRCh38) | `7674220` |
| ID | rsID или `.` если неизвестен | `rs28934574` |
| REF | Референсный аллель | `C` |
| ALT | Альтернативный аллель | `T` |
| QUAL | Фред-подобный скор качества варианта | `60.0` |
| FILTER | PASS или причина фильтрации | `PASS`, `LowQual` |
| INFO | Полуколонка-разделённые поля | `DP=45;AF=0.45;...` |

---

## INFO поля (в нашем демо VCF)

| Поле | Тип | Описание |
|------|-----|----------|
| `DP` | Integer | Суммарная глубина покрытия в позиции |
| `AF` | Float | Аллельная частота альтернативного аллеля в образце |
| `gnomAD_AF` | Float | Аллельная частота в gnomAD (популяционная) |
| `ClinVar` | String | Клиническое значение из ClinVar |

---

## VEP TSV — колонки

| Колонка | Описание |
|---------|----------|
| `Uploaded_variation` | ID варианта из VCF |
| `Location` | `CHROM:POS` или `CHROM:START-END` |
| `Allele` | Альтернативный аллель (VEP) |
| `Gene` | Ensembl Gene ID (ENSG...) |
| `Feature` | Ensembl Feature ID (ENST... для транскриптов) |
| `Feature_type` | Тип feature: Transcript, RegulatoryFeature, MotifFeature |
| `Consequence` | VEP consequence term(s), разделённые запятой |
| `cDNA_position` | Позиция в cDNA |
| `CDS_position` | Позиция в кодирующей последовательности |
| `Protein_position` | Позиция в белке |
| `Amino_acids` | Ref/Alt аминокислоты (для missense) |
| `Codons` | Ref/Alt кодоны |
| `Existing_variation` | rs ID из dbSNP и/или COSMIC |
| `Extra` | Дополнительные поля key=value через `;` |

### Extra — ключевые поля

| Ключ | Описание |
|------|----------|
| `IMPACT` | HIGH / MODERATE / LOW / MODIFIER |
| `SYMBOL` | Название гена (HGNC) |
| `CANONICAL` | YES если это canonical транскрипт |
| `BIOTYPE` | protein_coding, lncRNA, pseudogene и др. |
| `HGVSc` | Нотация в кДНК (c.XXXN>N) |
| `HGVSp` | Нотация в белке (p.AminoNNAamino) |
| `gnomADe_AF` | Аллельная частота в gnomAD экзомах |
| `ClinVar_CLNSIG` | ClinVar клиническое значение |
| `SIFT` | SIFT предсказание (tolerated / deleterious) |
| `PolyPhen` | PolyPhen-2 предсказание (benign / probably_damaging) |

---

## VEP Consequence Terms — иерархия

### HIGH impact (нарушают функцию белка)

| Consequence | Описание | Пример |
|-------------|----------|--------|
| `frameshift_variant` | Сдвиг рамки чтения (indel не кратный 3) | Ins/Del |
| `stop_gained` | Появление стоп-кодона (нонсенс) | p.Glu55* |
| `stop_lost` | Потеря стоп-кодона | p.Ter394Leu |
| `start_lost` | Потеря стартового метионина | p.Met1? |
| `splice_donor_variant` | Нарушение donor сайта сплайсинга (+1/+2) | c.934+1G>A |
| `splice_acceptor_variant` | Нарушение acceptor сайта сплайсинга (-1/-2) | c.935-2A>G |
| `transcript_ablation` | Делеция всего транскрипта | SV/CNV |

### MODERATE impact (изменяют белок, не разрушают)

| Consequence | Описание |
|-------------|----------|
| `missense_variant` | Замена аминокислоты (не stop) |
| `inframe_insertion` | Вставка аминокислоты (кратная 3) |
| `inframe_deletion` | Делеция аминокислоты (кратная 3) |
| `protein_altering_variant` | Изменение белка, не missense |

### LOW impact (минимальные функциональные последствия)

| Consequence | Описание |
|-------------|----------|
| `synonymous_variant` | Замена нуклеотида без изменения АК (молчащая) |
| `stop_retained_variant` | Замена в стоп-кодоне, стоп сохраняется |
| `splice_region_variant` | Изменение в 3–8 нуклеотидах от сайта сплайсинга |
| `start_retained_variant` | Изменение в стартовом кодоне, Met сохраняется |

### MODIFIER (влияние неизвестно или минимально)

| Consequence | Описание |
|-------------|----------|
| `intron_variant` | В интроне |
| `3_prime_UTR_variant` | В 3'UTR |
| `5_prime_UTR_variant` | В 5'UTR |
| `upstream_gene_variant` | До 5kb выше гена |
| `downstream_gene_variant` | До 5kb ниже гена |
| `intergenic_variant` | Между генами |
| `non_coding_transcript_variant` | В некодирующем транскрипте |

---

## SV VCF — дополнительные INFO поля

| Поле | Описание | Пример |
|------|----------|--------|
| `SVTYPE` | Тип структурного варианта | `DEL`, `DUP`, `INV`, `BND`, `INS` |
| `SVLEN` | Длина SV в bp (отрицательная для DEL) | `-20000` |
| `END` | Конечная позиция SV | `7694220` |
| `CIPOS` | Доверительный интервал для POS | `-10,10` |
| `CIEND` | Доверительный интервал для END | `-10,10` |
| `GENE` | Поражённые гены (нестандартное поле) | `TP53` |

### SVTYPE — расшифровка

| SVTYPE | Расшифровка | Клиническое значение |
|--------|-------------|----------------------|
| `DEL` | Deletion — потеря хромосомного материала | LOF, если поражает coding sequence |
| `DUP` | Duplication — усиление (tandem или inserted) | Может активировать онкогены |
| `INV` | Inversion — инверсия сегмента | Часто нарушает гены на границах |
| `BND` | Breakend — транслокация или сложная перестройка | Слияние генов (BCR-ABL, EML4-ALK) |
| `INS` | Insertion — вставка (в т.ч. мобильные элементы) | Может нарушать экзоны |

### Символьные ALT в SV VCF

| ALT | Описание |
|-----|----------|
| `<DEL>` | Делеция |
| `<DUP>` | Дупликация (танdemная) |
| `<INV>` | Инверсия |
| `<DUP:TANDEM>` | Тандемная дупликация |
| `<INS:ME>` | Вставка мобильного элемента |
| `N[chr:pos[` | BND: разрыв и соединение с chr:pos |

---

## Биологические базы данных

| База | Что хранит | URL |
|------|-----------|-----|
| **ClinVar** | Клинически интерпретированные варианты | ncbi.nlm.nih.gov/clinvar |
| **gnomAD** | Популяционные частоты (~730k экзомов) | gnomad.broadinstitute.org |
| **dbSNP** | Все известные полиморфизмы (rs ID) | ncbi.nlm.nih.gov/snp |
| **COSMIC** | Соматические мутации рака | cancer.sanger.ac.uk/cosmic |
| **OMIM** | Менделирующие заболевания | omim.org |
| **HGMD** | Мутации в наследственных болезнях | hgmd.cf.ac.uk |
