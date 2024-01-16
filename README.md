# sage-2.2.0
Next minor version of the Sage force field. **In development, do not use**

## Changes from Sage 2.1.0
All of these changes are a result of the [torsion multiplicity][tm] work
attempting to separate proper torsion parameters such that they apply only to
central bonds with a single multiplicity.

### Parameter changes
- t73 split into t73 and t73a
- t74 split into t74 and t74a
- t82 split into t82 and t82ax
- t83 split into t83 and t83ax
- t116 split into t116, t116a, t116b, and t116c
- t118 split into t118 and t118a
- t121 split into t121 and t121a
- t122 split into t122a, b, c, d, e, and f
- t130 split into t130a, b, c, and d
- t132 split into t132a, b, c, and d
- t133 split into t133 and t133a
- t142 split into t142a, b, c, d, e, and f
- t143 split into t143a, b, c, d, e, and f
- t157 split into t157 and t157a

### Data set changes
- Sage 2.1.0 optimization set augmented with [OpenFF multiplicity correction
  optimization set v1.0][tm-opt]
- Sage 2.1.0 torsion drive set augmented with [OpenFF multiplicity correction
  torsion drive data v1.1][tm-td]

### Environment changes
See [fb-196.yaml](fb-196.yaml) for the full fitting environment. Major changes
include:
- [ForceBalance][fb] 1.9.6, up from 1.9.3
- [openff-toolkit][offtk] 0.14.5, up from 0.10.6

[tm]: https://openforcefield.atlassian.net/wiki/spaces/FF/pages/2603909164/Torsion+multiplicity
[tm-opt]: https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-04-29-OpenFF-multiplicity-correction-optimization-set-v1.0
[tm-td]: https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-04-29-OpenFF-multiplicity-correction-torsion-drive-data
[fb]: https://github.com/leeping/forcebalance
[offtk]: https://github.com/openforcefield/openff-toolkit
