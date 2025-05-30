# Building the Random Effects Structure from $(r \mid f)$ in `lme4`

We aim to understand how `lme4` constructs the **relative covariance factor** $\mathbf{\Lambda}_\theta$ and the **random effects model matrix** $\mathbf{Z}$, given a model formula of the form $(r \mid f)$, where:

* $r$ is a **model formula** (e.g., `1` or `Time`), and
* $f$ is a **grouping factor** (e.g., `Subject`)



### Example: Random Intercepts and Slopes for Subject

Consider an example with $n = 10$ observations and a grouping factor `Subject` with $\ell = 3$ levels: A, B, and C. Suppose we include two random effects:

* A random intercept: $(1 \mid \text{Subject})$
* A random slope on `Time`: $(\text{Time} \mid \text{Subject})$



### Raw Random-Effects Model Matrices $X_i$

The expression $r$ is evaluated into an R model matrix $X_i$, of dimension $n \times p_i$, called the **raw random-effects model matrix** for the $i$th term.

* For $(1 \mid \text{Subject})$, we have:

$$
X_1 = \begin{bmatrix}
1 \\
1 \\
\vdots \\
1
\end{bmatrix}_{10 \times 1}
$$

* For $(\text{Time} \mid \text{Subject})$, we have:

$$
X_2 = \begin{bmatrix}
t_1 \\
t_2 \\
\vdots \\
t_{10}
\end{bmatrix}_{10 \times 1}
$$



### Grouping Factor Vector and Indicator Matrix

The expression $f$ evaluates to an R factor, called the **grouping factor**, which maps each observation to a group level. For the $i$th term, this is represented by a vector $i_i$ of **factor indices**, an $n$-vector with entries in $\{1, \dots, \ell_i\}$.

From this, we construct the **indicator matrix** $\mathbf{J}_i \in \mathbb{R}^{n \times \ell_i}$, whose columns correspond to group levels and whose rows indicate the group membership of each observation.

## Constructing the Model Matrix $Z$

Since the grouping factor `Subject` has three levels: A, B, and C. Then the indicator matrix $\mathbf{J}_i \in \mathbb{R}^{10 \times 3}$ is:

$$
\mathbf{J}_i = \begin{bmatrix}
1 & 0 & 0 \\
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \\
0 & 0 & 1 \\
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0
\end{bmatrix}
$$

Each row indicates the group (Subject A, B, or C) associated with that observation; each column corresponds to one group.

We now construct the matrix $\mathbf{Z}_i = \mathbf{J}_i \circ \mathbf{X}_i$, where $\circ$ denotes the **row-wise Kronecker product** (i.e., the **Khatri–Rao product**). This operation gives group-specific versions of each column of $\mathbf{X}_i$.

For the **random intercepts**, we have:

$$
Z_1 = \mathbf{J}_1 \circ X_1 \in \mathbb{R}^{10 \times 3}
$$

where

$$
Z_1 = \begin{bmatrix}
1 & 0 & 0 \\
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \\
0 & 0 & 1 \\
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0
\end{bmatrix}
$$

For the **random slopes on** `Time`, we have:

$$
Z_2 = \mathbf{J}_2 \circ X_2 \in \mathbb{R}^{10 \times 3}
$$

with

$$
Z_2 = \begin{bmatrix}
t_1 & 0 & 0 \\
t_2 & 0 & 0 \\
0 & t_3 & 0 \\
0 & t_4 & 0 \\
0 & 0 & t_5 \\
0 & 0 & t_6 \\
t_7 & 0 & 0 \\
0 & t_8 & 0 \\
0 & 0 & t_9 \\
t_{10} & 0 & 0
\end{bmatrix}
$$

The full **random-effects model matrix** $Z$ is then constructed by horizontally stacking the two blocks:

$$
Z = [Z_1 \mid Z_2] \in \mathbb{R}^{10 \times 6}
$$

$$
Z = \begin{bmatrix}
1 & 0 & 0 & t_1 & 0 & 0 \\
1 & 0 & 0 & t_2 & 0 & 0 \\
0 & 1 & 0 & 0 & t_3 & 0 \\
0 & 1 & 0 & 0 & t_4 & 0 \\
0 & 0 & 1 & 0 & 0 & t_5 \\
0 & 0 & 1 & 0 & 0 & t_6 \\
1 & 0 & 0 & t_7 & 0 & 0 \\
0 & 1 & 0 & 0 & t_8 & 0 \\
0 & 0 & 1 & 0 & 0 & t_9 \\
1 & 0 & 0 & t_{10} & 0 & 0
\end{bmatrix}
$$

Each group (Subject A, B, and C) has its own column for both the intercept and the slope.

## Constructing the Relative Covariance Factor $\mathbf{\Lambda}_\theta$

The matrix $\mathbf{\Lambda}_\theta$ encodes the **term-wise relative covariance factor** for the random effects. It is a **block-diagonal matrix** where each diagonal block corresponds to one random-effects term in the model. Each block, $\Lambda_i$, governs the covariance structure for all group levels in that term.


### Template Matrix $T_i$ and Parameter Vector $\theta$

Suppose for the random slope and intercept term $(1 + \text{Time} \mid \text{Subject})$, we have $p = 2$ random effects per group. The **template matrix** $T$ is a $2 \times 2$ lower-triangular matrix:

$$
T = \begin{bmatrix}
\theta_1 & 0 \\
\theta_2 & \theta_3
\end{bmatrix}
$$

This matrix encodes the random intercept variance ($\theta_1^2$), random slope variance ($\theta_3^2$), and intercept–slope covariance ($\theta_1 \theta_2$). In `lme4`, only the **lower-triangular elements** are parameterized to enforce symmetry and positive semi-definiteness.

The parameter vector $\theta$ contains these entries:

$$
\theta = (\theta_1, \theta_2, \theta_3)
$$

To ensure identifiability, `lme4` requires that the diagonal elements $\theta_1, \theta_3$ be **non-negative**.


### Block-Diagonal Structure of $\mathbf{\Lambda}_\theta$

For $\ell = 3$ groups (Subject A, B, C), each having the same $2 \times 2$ covariance structure, the full $\mathbf{\Lambda}_\theta$ is constructed by repeating the template matrix $T$ on the block diagonal:

$$
\mathbf{\Lambda}_\theta = \begin{bmatrix}
T & 0 & 0 \\
0 & T & 0 \\
0 & 0 & T
\end{bmatrix}
\in \mathbb{R}^{6 \times 6}
$$

Explicitly:

$$
\mathbf{\Lambda}_\theta = \begin{bmatrix}
\theta_1 & 0       & 0        & 0       & 0        & 0 \\
\theta_2 & \theta_3 & 0        & 0       & 0        & 0 \\
0        & 0       & \theta_1 & 0       & 0        & 0 \\
0        & 0       & \theta_2 & \theta_3 & 0        & 0 \\
0        & 0       & 0        & 0       & \theta_1 & 0 \\
0        & 0       & 0        & 0       & \theta_2 & \theta_3 \\
\end{bmatrix}
$$

This structure is **sparse** and **block-repetitive**.


### Covariance Matrix from $\mathbf{\Lambda}_\theta$

The covariance matrix for the full random effects vector $\mathbf{b}$ is:

$$
\Sigma_\theta = \sigma^2 \mathbf{\Lambda}_\theta \mathbf{\Lambda}_\theta^\top
$$

This implies that each group shares the same internal covariance structure $T T^\top$ for its associated random effects.

### Mapping the Parameter Vector $\theta$ into the Covariance Factor: Construction of `Lind`

The following is an analysis of how structured covariance matrices are implemented in the development branch [flexLambda](https://github.com/lme4/lme4/blob/flexLambda/R/reGenerators.R) of the `lme4` package.

The vector of random-effects covariance parameters, $\theta$, is compact: its length corresponds only to the number of *free* parameters in the model. For example, in the case of a diagonal covariance matrix with two random effects per group (e.g., an intercept and a slope), $\theta$ might contain just two elements: one for each standard deviation.

However, the relative covariance factor matrix $\Lambda_\theta$ — represented computationally in `lme4` as `Lambdat` — contains a full block-diagonal structure, with one lower-triangular block per level of the grouping factor. Each of these blocks is a Cholesky factor, and is typically sparse. Importantly, the matrix `Lambdat` is stored in **transposed form** and in **sparse column-major format**. This means only the non-zero elements in the **upper triangle** of each block are stored, and these entries are stacked column-by-column across all blocks.

In order to assign the appropriate values from $\theta$ into the correct entries of `Lambdat@x`, `lme4` uses an integer index vector named `Lind`. This vector has the same length as `Lambdat@x`, and serves as a mapping: the $i$-th element of `Lambdat@x` is set to the $j$-th element of $\theta$ if `Lind[i] = j`. That is,

$$
\texttt{Lambdat@x}[i] \leftarrow \theta[\texttt{Lind}[i]]
$$

This mapping ensures that shared parameters — such as a common variance or correlation term across all groups — are reused consistently.



#### Construction Phase

The construction of `Lambdat` and `Lind` begins with a **template Cholesky block** corresponding to a single level of the grouping factor. This block may be fully unstructured, diagonal, compound symmetric, or follow another pattern depending on the random-effects structure.

The Cholesky factor of this template is represented symbolically as a sparse matrix with placeholder entries (often created via `upper.tri()` or `chol()` on a covariance matrix). The vector of non-zero values in this template, say $\text{vec}(L_1)$, determines the **pattern of entries** in the full `Lambdat` matrix. The number of parameters required for a single block is:

$$
m = \frac{p(p + 1)}{2}
$$

where $p$ is the number of random effects per group. These parameters are indexed from $1$ to $m$ and stored in a vector `Lind_template` of length equal to the number of non-zero entries in $\text{vec}(L_1)$, typically ordered in column-major (Fortran-style) fashion.

Once the pattern for a single group is known, it is replicated for each level of the grouping factor:

```r
Lind <- rep(Lind_template, times = ell)
````

where $\ell$ is the number of grouping levels. The matrix `Lambdat` itself is constructed using:

```r
Lambdat <- bdiag(replicate(ell, list(<template sparse block>)))
```

The resulting `Lambdat` matrix has the correct sparsity pattern, but the actual values in its non-zero entries are initially set to 1s or some placeholder values. These will be updated at runtime.



#### Update Phase

During model fitting, the optimizer proposes new values for the parameter vector $\theta$. To reflect these new values in the model, the corresponding covariance factor must be updated.

This is achieved via the function `updateLambdatx(theta)`, which simply returns:

```r
theta[Lind]
```

This vector is then assigned directly into the sparse matrix object:

```r
Lambdat@x <- updateLambdatx(theta)
```

At this point, `Lambdat@x` holds the actual numeric values for the Cholesky factor of the relative covariance matrix — repeated block-wise across all grouping levels. The structure of `Lambdat` remains fixed, while its contents are dynamically updated with each new proposal of \$\theta\$ during optimization.



As a result, `Lind` serves as the glue between the compact parameter vector $\theta$ and the expanded block-sparse matrix `Lambdat`. 


### Example: Diagonal Covariance Structure

To illustrate how the `theta` vector and the `Lind` index map into the transposed relative covariance factor `Lambdat`, consider the simple model we used above with a diagonal covariance structure and two random effects per group, a random intercept and a random slope on `Days`:

```r
Reaction ~ Days + (Days | Subject)
````

Assume there are $\ell = 3$ levels of the grouping factor `Subject`, and $p = 2$ random effects per group. Under a diagonal structure (i.e., no correlation between intercept and slope), each group-level Cholesky block is a $2 \times 2$ diagonal matrix:

$$
L_i = \begin{bmatrix}
\sigma_1 & 0 \\
0 & \sigma_2
\end{bmatrix}
$$

The parameter vector is:

```r
theta <- c(1.5, 0.8)  # sigma_1 = 1.5 (intercept), sigma_2 = 0.8 (slope)
```

This vector contains the standard deviations of the intercept and slope effects, and its length is $m = p = 2$.



#### Construction of `Lind`

To populate `Lambdat@x` with the correct entries from `theta`, we define an index vector `Lind` of length $2 \times 3 = 6$ (i.e., 2 values per group, 3 groups):

```r
Lind <- rep(1:2, times = 3) # c(1, 2, 1, 2, 1, 2)
```

This tells `lme4` to apply`theta[1]`($\sigma_1$) to all intercept terms and `theta[2]`($\sigma_2$) to all slope terms

---

#### Construction of `Lambdat`

We now construct a symbolic `Lambdat` matrix with the appropriate sparsity pattern:

```r
Lambdat <- bdiag(replicate(3, Diagonal(x = c(1, 1))))
```

This gives a $6 \times 6$ block-diagonal matrix with $2 \times 2$ identity matrices on the diagonal (one per group), initially filled with placeholders.

The actual numerical entries of the matrix are set during the update phase:

```r
Lambdat@x <- theta[Lind]
# c(1.5, 0.8, 1.5, 0.8, 1.5, 0.8)
```

The resulting transposed covariance factor matrix is:

$$
\Lambda_\theta^\top = 
\begin{bmatrix}
1.5 & 0   & 0   & 0   & 0   & 0   \\
0   & 0.8 & 0   & 0   & 0   & 0   \\
0   & 0   & 1.5 & 0   & 0   & 0   \\
0   & 0   & 0   & 0.8 & 0   & 0   \\
0   & 0   & 0   & 0   & 1.5 & 0   \\
0   & 0   & 0   & 0   & 0   & 0.8
\end{bmatrix}
$$

Here, each $2 \times 2$ diagonal block corresponds to one level of the grouping factor (`Subject`), with standard deviations applied according to their position in the `theta` vector.




