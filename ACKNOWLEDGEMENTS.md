# Acknowledgements and provenance / 致谢与来源

## 中文

`ggideogram` 是一个独立的 ggplot2-first 重构项目，但它建立在
[RIdeogram](https://github.com/TickingClock1992/RIdeogram) 的原创工作之上。
没有 RIdeogram 对全基因组 ideogram 可视化的设计、实现和开放发布，就不会有本项目。

我们郑重感谢 RIdeogram 作者：Zhaodong Hao、Dekang Lv、Ying Ge、Jisen Shi、
Dolf Weijers、Guangchuang Yu 和 Jinhui Chen。

本项目从 RIdeogram 继承或借鉴了：

- 染色体胶囊、圆帽、着丝粒及全基因组排列的几何思想；
- `human_karyotype`、`liriodendron_karyotype`、`gene_density`、
  `LTR_density` 和 `Random_RNAs_500` 示例数据；
- `GFFex()`、基因组窗口划分及经典示例的行为参考；
- 使用 ideogram 映射全基因组标记、密度和区间信息的工作范式。

`ggideogram` 随后的布局、坐标变换、原生 ggplot2 marker、声明式轨道、完整图形 inset、
原生染色体 x/y guide、响应式物理尺寸和组合协议均为独立重构。本仓库不是 RIdeogram
作者维护的官方后继版本。

使用本包或上述继承数据时，请同时引用 `ggideogram` 和 RIdeogram 原始论文：

> Hao Z, Lv D, Ge Y, Shi J, Weijers D, Yu G, Chen J. (2020).
> RIdeogram: drawing SVG graphics to visualize and map genome-wide data on
> the idiograms. *PeerJ Computer Science*, 6:e251.
> https://doi.org/10.7717/peerj-cs.251

RIdeogram 可从 [GitHub](https://github.com/TickingClock1992/RIdeogram) 和
[CRAN](https://cran.r-project.org/package=RIdeogram) 获取。原包与本项目均采用
Artistic License 2.0。

## English

`ggideogram` is an independent ggplot2-first refactor built on the original
work of [RIdeogram](https://github.com/TickingClock1992/RIdeogram). This project
would not exist without RIdeogram's design, implementation, and open release
for genome-wide ideogram visualization.

We gratefully acknowledge the RIdeogram authors: Zhaodong Hao, Dekang Lv,
Ying Ge, Jisen Shi, Dolf Weijers, Guangchuang Yu, and Jinhui Chen.

The following originated in, or were informed by, RIdeogram:

- the chromosome capsule, rounded-cap, centromere, and whole-genome layout
  concepts;
- the `human_karyotype`, `liriodendron_karyotype`, `gene_density`,
  `LTR_density`, and `Random_RNAs_500` example datasets;
- behavioural references for `GFFex()`, genomic windows, and classic examples;
- the workflow of mapping genome-wide markers, densities, and intervals onto
  idiograms.

The dimensionless layout, coordinate transforms, native ggplot2 markers,
declarative tracks, complete-plot insets, native chromosome x/y guides,
physical-size scaling, and composition protocols in `ggideogram` are an
independent redesign. This repository is not an official successor maintained
by the RIdeogram authors.

When using this package or the inherited datasets above, please cite both
`ggideogram` and the original RIdeogram paper:

> Hao Z, Lv D, Ge Y, Shi J, Weijers D, Yu G, Chen J. (2020).
> RIdeogram: drawing SVG graphics to visualize and map genome-wide data on
> the idiograms. *PeerJ Computer Science*, 6:e251.
> https://doi.org/10.7717/peerj-cs.251

RIdeogram is available from
[GitHub](https://github.com/TickingClock1992/RIdeogram) and
[CRAN](https://cran.r-project.org/package=RIdeogram). Both the upstream package
and this project are released under the Artistic License 2.0.
