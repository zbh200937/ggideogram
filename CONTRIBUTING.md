# Contributing / 贡献指南

## 中文

欢迎 issue 和 pull request。提交前请：

1. 先搜索现有 issue，说明问题对应的染色体方向、轨道类型和输出设备；
2. 提供最小可运行示例，并优先使用包内小型数据；
3. 不上传原始测序数据、个人目录、令牌、缓存或大型中间文件；
4. 为行为变化补充紧邻源码的 `testthat` 测试；
5. 运行：

```sh
Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat")'
R CMD build . --no-build-vignettes
R CMD check --no-manual ggideogram_*.tar.gz
```

涉及布局或组合协议的变更还应先阅读
[ARCHITECTURE.md](ARCHITECTURE.md)。请保留 RIdeogram 的来源说明和数据署名。

## English

Issues and pull requests are welcome. Before submitting:

1. search existing issues and identify the chromosome orientation, track type,
   and output device involved;
2. provide a minimal reproducible example, preferably with bundled data;
3. never upload raw sequencing data, personal paths, tokens, caches, or large
   intermediate files;
4. add focused `testthat` coverage for behavioural changes;
5. run:

```sh
Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat")'
R CMD build . --no-build-vignettes
R CMD check --no-manual ggideogram_*.tar.gz
```

Changes to layout or composition contracts should also follow
[ARCHITECTURE.md](ARCHITECTURE.md). Preserve the RIdeogram provenance and data
attribution.
