# SH + ILO
- Fast and accurate hinge detection algorithm by comparing two different conformationf of the same protein

## Preliminaries
### RMSD [1]
- $n$: the number of aligned residues
- $A = \{\vec{a_1}, \vec{a_2},\dots \vec{a_n}\}$
- $B = \{\vec{b_1}, \vec{b_2},\dots \vec{b_n}\}$
- $\mathcal{R}$: a set of possible rotation matrix
- $\mathcal{V}$: a set of possible translation vectors
- $A[i .. j]$: Substructure of $A$ for $i, j \in \mathbb{N}_1^n$ such that $i \leq j$.
- $sd(A, B) = \min_{\vec{v} \in \mathcal{V}, R \in \mathcal{R}}||\vec{a_i} - R(\vec{b} - \vec{v})||^2$

$$
RMSD(A, B) = \min_{R \in \mathcal{R}, \vec{v} \in \mathcal{V}}\sqrt{\frac{1}{n}\sum_{i=1}^n |\vec{a_i} - R(\vec{b_i} - \vec{v})|}
$$
### $RMSDh^{(k)}$ [2]

$$
RMSDh^{(k)}(A, B) = \min_{1 \leq h_1 \leq h_2 \leq \dots \leq h_k \leq n} \sqrt{\frac{1}{n}\sum_{i=1}^{k+1}sd(A[h_{i-1}..h_i - 1], B[h_i.. h_i-1])}
$$

#### Row-minimum index problem
Let $\mathbb{R}$ denote a set of real numbers.
Let $\mathbb{N}$ denote a set of natural numbers.
For a natural number $n$, let $\mathbb{N}_1^n = \{x \in \mathbb{N} | 1 \leq x \leq n\}$.
Let $I_f(i) = \argmin_{j \in \mathbb{N}_1^n}f(i, j)$.
- Given a function $f : \mathbb{N}_1^n \times \mathbb{N}_1^n \rightarrow \mathbb{R}$, compute $I_f(i)$ for all integer $i \in \mathbb{N}_1^n$.

#### Dynamic Programming (DP)

The formulation of DP to solve the row-minimum indexed problem is as follows.

$$
D(i, r, A, B) = \min_{1 \leq h_1 \leq h_2 \leq \dots \leq h_r < i}\sum_{l=1}^{r+1}sd(A[h_{l-1}..h_l-1], B[h_{l-1}..h_l-1])
$$

$$
D(i, r, A, B) = \min_{r \leq j \leq i}\{D(j, r-1, A, B) + sd(A[j..i], B[j..i])\}
$$
The time complexity for solving the above problem is $O(kn^2)$.
By using SMAWK algorithm [3], the above problem can be solved in $O(kn)$.

### Experiments

Python analysis scripts require Python 3.11 or 3.12 and Poetry 2.2 or later.
Install the locked dependencies with `poetry install --no-root`.
The Python minimum is required by the security-fixed Pillow 13 dependency.

You can enter the virtual environment by the following commands
```bash
docker compose build --pull
docker compose up -d
docker container exec -it shilo bash
```

You can see the experimental result in the rmsdh_result directory.
```bash
cd src
unzip rmsdh_result.zip
```

You can reproduce the experimental results by the following code
```bash
cd src
unzip coord_csv.zip
bash speed_comparison_rmsdhk.sh
bash ablation_study.sh
```

You can also compute $RMSDh^{(k)}$ of two protein pairs by the following command
```bash
g++ execute_sh_ilo.cpp -o execute_sh_ilo -std=c++14 -O3
./execute_sh_ilo pdb3hvp.pdb pdb4hvp.pdb A A 2
```

Both command-line programs require a positive hinge count smaller than the
number of aligned CA residues. Missing files, malformed ATOM records, empty
chains, and invalid hinge counts exit with an error before computation.

### Security checks

Audit every package in `poetry.lock`, including dependencies for other Python
versions and platforms, using [pip-audit](https://github.com/pypa/pip-audit):

```bash
uv run --no-project --python 3.12 --with pip-audit python scripts/audit_dependencies.py
```

The container uses [Ubuntu 24.04 LTS](https://ubuntu.com/about/release-cycle).
Rebuild with `docker compose build --pull --no-cache` to pick up new OS fixes.
Run the C++ input regression tests with AddressSanitizer and
UndefinedBehaviorSanitizer from the repository root:

```bash
docker build --pull -t shilo-security-check .
docker run --rm --network none \
  --mount "type=bind,source=$PWD,target=/workspace,readonly" \
  --workdir /workspace shilo-security-check bash scripts/test_security.sh
```


### References
- [1] Kabsch, Wolfgang (1978). "A discussion of the solution for the best rotation to relate two sets of vectors". Acta Crystallographica. A34 (5): 827–828.
- [2] Shibuya, Tetsuo. "Fast hinge detection algorithms for flexible protein structures." IEEE/ACM Transactions on Computational Biology and Bioinformatics 7.2 (2008): 333-341.
- [3] Aggarwal, Alok, et al. "Geometric applications of a matrix searching algorithm." Proceedings of the second annual symposium on Computational geometry. 1986.
