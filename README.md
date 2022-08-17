## T-Cell Receptor-Peptide Interaction Prediction with Physical Model Augmented Pseudo-Labeling

The **[paper](https://dl.acm.org/doi/10.1145/3534678.3539075)** by [Yiren Jian](https://cs.dartmouth.edu/~yirenjian/), [Erik Kruus](https://www.nec-labs.com/research/machine-learning/people/erik-kruus/) and [Martin Renqiang Min](https://www.cs.toronto.edu/~cuty/) is accepted to ACM SIGKDD (KDD) 2022. The work is done while Yiren doing an internship at **[NEC Labs America](https://www.nec-labs.com/)**.

- [ ] paper
- [x] code
- [x] data

```bibtex
@inproceedings{jian-etal-2022-Tcell,
  author = {Jian, Yiren and Kruus, Erik and Min, Martin Renqiang},
  title = {T-Cell Receptor-Peptide Interaction Prediction with Physical Model Augmented Pseudo-Labeling},
  year = {2022},
  isbn = {9781450393850},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3534678.3539075},
  doi = {10.1145/3534678.3539075},
  booktitle = {Proceedings of the 28th ACM SIGKDD Conference on Knowledge Discovery and Data Mining},
  pages = {3090–3097},
  numpages = {8},
  keywords = {pseudo-labeling, T-cell receptors, docking energy, peptide recognition, deep neural network, physical modeling},
  location = {Washington DC, USA},
  series = {KDD '22}
}
```


### Requirements and Related Works

Please see `PhyAugmentation/README.md` for installation instructions.

Our work is built on several prior works including: [ERGO](https://github.com/IdoSpringer/ERGO) and [HDOCK](http://hdock.phys.hust.edu.cn/).

```bibtex
@article{springer2020prediction,
  title={Prediction of specific TCR-peptide binding from large dictionaries of TCR-peptide pairs},
  author={Springer, Ido and Besser, Hanan and Tickotsky-Moskovitz, Nili and Dvorkin, Shirit and Louzoun, Yoram},
  journal={Frontiers in immunology},
  pages={1803},
  year={2020},
  publisher={Frontiers}
}
```

```bibtex
@article{yan2020hdock,
  title={The HDOCK server for integrated protein--protein docking},
  author={Yan, Yumeng and Tao, Huanyu and He, Jiahua and Huang, Sheng-You},
  journal={Nature protocols},
  volume={15},
  number={5},
  pages={1829--1852},
  year={2020},
  publisher={Nature Publishing Group}
}
```
