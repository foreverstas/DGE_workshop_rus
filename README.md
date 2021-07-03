**Это вольный перевод семинара [DGE_workshop](https://github.com/hbctraining/DGE_workshop). В процессе перевода никаких изменений в содержание курса не вводилось. К ключевым терминам прикреплены ссылки на интернет ресурсы с пояснениями**

# Differential gene expression workshop

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | [Introduction to R](https://hbctraining.github.io/Intro-to-R/) | 1.5-day workshop (~10 hours of trainer-led time)|

### Описание

Репозиторий содержит обучающие материалы **1.5-дневного** практического семинара **Введение в анализ дифференциальной экспрессии генов, DGE**. Семинар познакомит участников с анализом дифференциальной экспрессии генов по количественным результатам RNA-seq с использованием [R](https://starnew.inp.nsk.su/~baldin/DataAnalysis/R/R-01-intro.pdf)/[RStudio](https://ru.wikipedia.org/wiki/RStudio). Семинар проведет участников через выполнение рабочего процесса анализа дифференциальной экспрессии генов на данных подсчета РНК-секвенирований с использованием R / RStudio. Требуются рабочие знания языка программирования R или завершение [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/).

### Цели обучения

- Контроль качества (QC) численных данных с использованием Метода Главных Компонент [(PCA)](https://www.youtube.com/watch?v=_UVHneBUBW0) и [иерархической кластеризации](https://www.youtube.com/watch?v=_UVHneBUBW0)
- Использование DESeq2 для определения списка достоверно различающихися генов
- Визуализация паттернов экспрессии [дифференциально экспрессирующихся генов](https://ru.wikipedia.org/wiki/%D0%9A%D0%BE%D0%BB%D0%B8%D1%87%D0%B5%D1%81%D1%82%D0%B2%D0%B5%D0%BD%D0%BD%D1%8B%D0%B9_%D0%B0%D0%BD%D0%B0%D0%BB%D0%B8%D0%B7_%D1%8D%D0%BA%D1%81%D0%BF%D1%80%D0%B5%D1%81%D1%81%D0%B8%D0%B8_%D0%B3%D0%B5%D0%BD%D0%BE%D0%B2)
- Функциаланьный анализ списков генов с использованием инструментов языка программирования R

> These materials are developed for a trainer-led workshop, but also amenable to self-guided learning.

### Уроки

Ссылки на лекции:

* [Click here for schedule using Salmon count matrix](./lessons/01_DGE_setup_and_overview.html)
* [Click here for schedule using FeatureCounts count matrix](schedule/1.5-day.md)


### Требуемые программы

1. Загрузите наиболее свежие версии R и RStudio на устройство [Инструкция на YouTube](https://www.youtube.com/watch?v=xct_zaU5zL8):

 - [R](https://cran.r-project.org/)
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

2. Установите следующие пакеты как показано в инструкции ниже.

> **NOTE:**  При установке пакетов, при необходимости сделать выбор (a/s/n) или (y/n), вводите “a” или "y", но знайте, что в таком случае на установку потребуется больше времени.

(a) Установите на устройство **пакеты** из репозитория **CRAN**. Вам не нужно обращаться к странице CRAN в интернете, просто используйте эти функции для установки каждого пакета по очереди:

```r
install.packages("имя_первого_пакета_в_кавычках")
install.packages("имя_второго_пакета_в_кавычках")
и т.д. ...
```

Пакеты для установки из CRAN (учтите, что названия пакетов чувствительны к регистру!):

* BiocManager
* RColorBrewer
* pheatmap
* ggrepel
* devtools
* tidyverse


(b) Установите указанные ниже **пакеты** из репозитория **Bioconductor**, используя функцию `BiocManager::install()` 7 раз для 7 пакетов:

```r
BiocManager::install("имя_первого_пакета_в_кавычках")
BiocManager::install("имя_второго_пакета_в_кавычках")

```

Имена пакетов из Bioconductor (все имена чувствительны к регистру!):

* DESeq2
* clusterProfiler
* DOSE
* org.Hs.eg.db
* pathview
* DEGreport
* EnsDb.Hsapiens.v86
* AnnotationHub
* ensembldb


3. В конце проверьте, что все пакеты установлены успешно, просто загрузив их друг за другом с помощью функции library().  

```r
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)
```

4. Как только загрузили пакеты, запустите функцию sessionInfo().  

```r
sessionInfo()
```



****

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
