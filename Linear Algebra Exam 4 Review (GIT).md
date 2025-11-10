---
Title: "Linear Algebra Exam 4 Review"
layout: default
permalink: /
---

<!-- When publishing to GitHub Pages, rename this file to `index.md` and place it at the repo root (for username.github.io) or in `/docs` (for any repo). -->

# Linear Algebra Exam 4 Review

> [Home](/) · [Download as Markdown](./Linear%20Algebra%20Exam%204%20Review%20(GIT).md)

---

<span style="color:#8B008B; font-weight:bold;">## Part 1: Gram–Schmidt Process (QR Factorization) ✨</span>

The **Gram–Schmidt process** is a method for orthogonalizing a set of vectors in an inner product space, which leads to the QR factorization of a matrix. Given a set of linearly independent vectors \(\{ \mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n \}\), the goal is to produce an orthonormal set \(\{ \mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n \}\) such that each \(\mathbf{u}_i\) is orthogonal to all previous \(\mathbf{u}_j\) for \(j < i\), and each has unit length.

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
This means we start with vectors that might not be perpendicular or have length 1, and we transform them into a set of vectors that are all perpendicular to each other (orthogonal) and have length 1 (normalized). This is useful for many applications like simplifying matrix computations.

---

### The Process:

1. Set \(\mathbf{u}_1 = \mathbf{v}_1\).

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
We begin by taking the first vector as is because there are no previous vectors to make it orthogonal to. Later, we will normalize it to have length 1.

2. For each subsequent vector \(\mathbf{v}_k\), subtract the projections onto all previous \(\mathbf{u}_j\):

\[
\mathbf{w}_k = \mathbf{v}_k - \sum_{j=1}^{k-1} \text{proj}_{\mathbf{u}_j}(\mathbf{v}_k),
\]

where the projection is:

\[
\text{proj}_{\mathbf{u}_j}(\mathbf{v}_k) = \frac{\mathbf{v}_k \cdot \mathbf{u}_j}{\mathbf{u}_j \cdot \mathbf{u}_j} \mathbf{u}_j.
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
Here, we remove from \(\mathbf{v}_k\) the parts that lie in the directions of the previous \(\mathbf{u}_j\) vectors. This subtraction ensures that the new vector \(\mathbf{w}_k\) is orthogonal (perpendicular) to all the previous \(\mathbf{u}_j\). The projection tells us how much of \(\mathbf{v}_k\) points in the direction of \(\mathbf{u}_j\), so subtracting it removes that component.

3. Normalize \(\mathbf{w}_k\) to get \(\mathbf{u}_k\):

\[
\mathbf{u}_k = \frac{\mathbf{w}_k}{\|\mathbf{w}_k\|}.
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
After making \(\mathbf{w}_k\) orthogonal to the previous vectors, we scale it to have length 1. This is important because having unit vectors simplifies calculations and forms the orthonormal set.

---

### Example (matches board):

Given the **3×2 matrix**
\[
A=\begin{bmatrix}
1 & 1\\
1 & 0\\
1 & -1
\end{bmatrix}=\begin{bmatrix}a_1 & a_2\end{bmatrix},
\quad a_1=\begin{bmatrix}1\\1\\1\end{bmatrix},\ 
a_2=\begin{bmatrix}1\\0\\-1\end{bmatrix}.
\]

---

**Step 1 — Normalize \(a_1\):**
\[
\|a_1\|=\sqrt{1^2+1^2+1^2}=\sqrt{3},\qquad
q_1=\dfrac{a_1}{\|a_1\|}=\dfrac{1}{\sqrt{3}}\!\begin{bmatrix}1\\1\\1\end{bmatrix}.
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
We find the length (or magnitude) of \(a_1\) by taking the square root of the sum of the squares of its components. This length tells us how long the vector is. Then, we divide each component by this length to create a unit vector \(q_1\) pointing in the same direction as \(a_1\), but with length 1.

---

**Step 2 — Orthogonalize \(a_2\) against \(q_1\):**
\[
\operatorname{proj}_{q_1}(a_2)=(q_1^{\mathsf T}a_2)\,q_1
=\left(\tfrac{1}{\sqrt{3}}[1\ 1\ 1]\!\begin{bmatrix}1\\0\\-1\end{bmatrix}\right)\!q_1
=\left(\tfrac{1-0-1}{\sqrt{3}}\right)q_1=0.
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
We compute the projection of \(a_2\) onto \(q_1\) to find the part of \(a_2\) that points in the direction of \(q_1\). The dot product \(q_1^{\mathsf T}a_2\) measures how much \(a_2\) aligns with \(q_1\). Since this projection is zero, it means \(a_2\) is already orthogonal (perpendicular) to \(q_1\), so no subtraction is needed.

Hence \(u_2=a_2\).

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
Because the projection is zero, \(a_2\) is already orthogonal to \(q_1\), so we keep \(a_2\) as is for the next step.

---

**Step 3 — Normalize to get \(q_2\):**
\[
\|u_2\|=\sqrt{1^2+0^2+(-1)^2}=\sqrt{2},\qquad
q_2=\dfrac{u_2}{\|u_2\|}=\dfrac{1}{\sqrt{2}}\!\begin{bmatrix}1\\0\\-1\end{bmatrix}.
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
We find the length of \(u_2\) and then divide each component by this length to create a unit vector \(q_2\). This ensures \(q_2\) has length 1 and is orthogonal to \(q_1\).

---

**Step 4 — Assemble \(Q\) and compute \(R=Q^{\mathsf T}A\):**
\[
Q=\begin{bmatrix}q_1 & q_2\end{bmatrix}
=\begin{bmatrix}
\tfrac{1}{\sqrt{3}} & \tfrac{1}{\sqrt{2}}\\[2pt]
\tfrac{1}{\sqrt{3}} & 0\\[2pt]
\tfrac{1}{\sqrt{3}} & -\tfrac{1}{\sqrt{2}}
\end{bmatrix},
\]
\[
R=\begin{bmatrix}
r_{11} & r_{12}\\
0 & r_{22}
\end{bmatrix}
=\begin{bmatrix}
q_1^{\mathsf T}a_1 & q_1^{\mathsf T}a_2\\
0 & q_2^{\mathsf T}a_2
\end{bmatrix}
=\begin{bmatrix}
\sqrt{3} & 0\\
0 & \sqrt{2}
\end{bmatrix}.
\]

### Compute entries of \(R\) (numbers) — and what \(^{\mathsf T}\) means

\[
r_{11} = q_1^{\mathsf T} a_1 = \frac{1}{\sqrt{3}}[1\ 1\ 1]\!\begin{bmatrix}1\\1\\1\end{bmatrix}
= \frac{1+1+1}{\sqrt{3}} = \frac{3}{\sqrt{3}} = \sqrt{3}.
\]
\[
r_{12} = q_1^{\mathsf T} a_2 = \frac{1}{\sqrt{3}}[1\ 1\ 1]\!\begin{bmatrix}1\\0\\-1\end{bmatrix}
= \frac{1+0-1}{\sqrt{3}} = 0.
\]
\[
r_{22} = q_2^{\mathsf T} a_2 = \frac{1}{\sqrt{2}}[1\ 0\ -1]\!\begin{bmatrix}1\\0\\-1\end{bmatrix}
= \frac{1+0+1}{\sqrt{2}} = \frac{2}{\sqrt{2}} = \sqrt{2}.
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span> The little superscript \(^{\mathsf T}\) means **transpose**: it turns a column into a row so you can take a **dot product**. Each \(r_{ij}\) is the amount of \(a_j\) in the direction of \(q_i\) (a projection).

---

<span style="color:#800080; font-weight:bold;">💜 Conclusion:</span> \(\boxed{A=QR}\) with the \(Q\) and \(R\) above (orthonormal columns, upper–triangular \(R\)).

---

### Summary:
- Columns of \(Q\) are \(q_1=\frac{1}{\sqrt{3}}[1,1,1]^{\mathsf T}\) and \(q_2=\frac{1}{\sqrt{2}}[1,0,-1]^{\mathsf T}\).
- \(R=\begin{bmatrix}\sqrt{3}&0\\0&\sqrt{2}\end{bmatrix}\).
- Matches the board calculation; no projection term for \(a_2\) because \(q_1^{\mathsf T}a_2=0\).

---

---

<span style="color:#8B008B; font-weight:bold;">## Part 2: Eigenvalues and Diagonalization 📐</span>

---

### Example 1: Characteristic Polynomial Derivation (3×3 Symmetric Matrix)

This example demonstrates the step-by-step derivation of the characteristic polynomial for a different matrix from lecture.

---

Putting it all together:

\[
6 - 11 \lambda + 6 \lambda^2 - \lambda^3 - \color{red}{(2-\lambda - 1)} - \color{red}{(1 - \lambda)} = 0,
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
This equation comes from expanding the determinant of \((A - \lambda I)\), which is used to find eigenvalues. The goal is to find values of \(\lambda\) that make this determinant zero.

Simplify the subtracted terms:

\[
\color{red}{(2-\lambda - 1)} = \color{green}{1 - \lambda}, \quad \text{so subtracting gives } \color{red}{- (1 - \lambda)} = \color{green}{-1 + \lambda}.
\]

Similarly,

\[
\color{red}{- (1 - \lambda)} = \color{green}{-1 + \lambda}.
\]

So total subtraction is:

\[
\color{green}{- (1 - \lambda) - (1 - \lambda) = -2 + 2 \lambda.}
\]

---

Therefore, characteristic polynomial:

\[
6 - 11 \lambda + 6 \lambda^2 - \lambda^3 \color{green}{- 2 + 2 \lambda} = 0,
\]

which simplifies to:

\[
(6 - 2) + (-11 \lambda + 2 \lambda) + 6 \lambda^2 - \lambda^3 = 0,
\]

\[
\color{blue}{4 - 9 \lambda + 6 \lambda^2 - \lambda^3 = 0.}
\]

Rewrite:

\[
\color{blue}{-\lambda^3 + 6 \lambda^2 - 9 \lambda + 4 = 0.}
\]

Multiply both sides by \(-1\):

\[
\color{blue}{\lambda^3 - 6 \lambda^2 + 9 \lambda - 4 = 0.}
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
This polynomial equation is called the characteristic polynomial. The roots (solutions) of this polynomial are the eigenvalues of the matrix. Finding these roots tells us the special values \(\lambda\) where the matrix behaves in a particular way.

---

### Example 2 (Board Problem, corrected):  \( A = \begin{bmatrix}4&1&0\\0&4&0\\0&0&2\end{bmatrix} \)

---

### 1. Characteristic Polynomial

\[
\det(A - \lambda I) = \det \begin{bmatrix}
4 - \lambda & 1 & 0 \\
0 & 4 - \lambda & 0 \\
0 & 0 & 2 - \lambda
\end{bmatrix}
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
We subtract \(\lambda I\) from \(A\) because we want to find values of \(\lambda\) for which \((A - \lambda I)\) is singular (i.e., its determinant is zero). This condition defines the eigenvalues.

Since the matrix is upper-triangular, the determinant is the product of diagonal entries:

\[
= (4 - \lambda)(4 - \lambda)(2 - \lambda) = (4 - \lambda)^2 (2 - \lambda).
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
For upper-triangular matrices, the determinant is the product of the diagonal elements. This simplifies the calculation significantly.

---

### 2. Eigenvalues

Set the characteristic polynomial equal to zero:

\[
(4 - \lambda)^2 (2 - \lambda) = 0.
\]

So eigenvalues are:

\[
\boxed{
\lambda_1 = 4 \quad (\text{multiplicity } 2), \quad \lambda_2 = 2 \quad (\text{multiplicity } 1).
}
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
The eigenvalues are the roots of the characteristic polynomial. Here, \(\lambda=4\) is a root twice (multiplicity 2), and \(\lambda=2\) is a root once (multiplicity 1). Multiplicity indicates how many times an eigenvalue repeats as a root.

---

### 3. Eigenvectors

---

#### For \(\lambda = 4\):

Solve:

\[
(A - 4I) \mathbf{v} = 0,
\]

\[
\begin{bmatrix}
0 & 1 & 0 \\
0 & 0 & 0 \\
0 & 0 & -2
\end{bmatrix}
\begin{bmatrix} x \\ y \\ z \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}.
\]

This gives the system:

\[
0 \cdot x + 1 \cdot y + 0 \cdot z = 0 \implies y = 0,
\]
\[
0 = 0 \quad (\text{no information}),
\]
\[
0 \cdot x + 0 \cdot y - 2 z = 0 \implies z = 0.
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
Plugging in \(\lambda=4\) into \((A - \lambda I)\) and solving \((A - 4I)\mathbf{v} = 0\) finds vectors \(\mathbf{v}\) that, when multiplied by \(A\), only get scaled by 4. This means these vectors don't change direction under \(A\), only length.

---

Eigenvectors satisfy:

\[
y = 0, \quad z = 0,
\]

and \(x\) is free. So eigenvectors are:

\[
\mathbf{v} = \begin{bmatrix} x \\ 0 \\ 0 \end{bmatrix} = x \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}.
\]

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
Since \(x\) can be any number, the eigenvectors for \(\lambda=4\) form a one-dimensional space along the \(x\)-axis.

---

#### For \(\lambda = 2\):

Solve:

\[
(A - 2I) \mathbf{v} = 0,
\]

\[
\begin{bmatrix}
2 & 1 & 0 \\
0 & 2 & 0 \\
0 & 0 & 0
\end{bmatrix}
\begin{bmatrix} x \\ y \\ z \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}.
\]

The system is:

\[
2x + y = 0,
\]
\[
2y = 0 \implies y = 0,
\]
\[
0 = 0 \quad (\text{no information}).
\]

---

From second equation: \( y = 0 \).

From first equation:

\[
2x + 0 = 0 \implies x = 0.
\]

Third coordinate \(z\) is free.

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
Here, plugging in \(\lambda=2\) finds vectors that only get scaled by 2 under \(A\). The system restricts \(x\) and \(y\) to zero, but \(z\) can be any value, meaning eigenvectors lie along the \(z\)-axis.

---

Eigenvectors:

\[
\mathbf{v} = \begin{bmatrix} 0 \\ 0 \\ z \end{bmatrix} = z \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}.
\]

<span style="color:#DC143C; font-weight:bold;">⚠️ Professor slip corrected:</span> \(E_{2}\) is 1‑dimensional (only \([0,0,1]^{\mathsf T}\)); it is **not** 2‑dimensional.

---

### 4. Multiplicities and Eigenspaces

| Eigenvalue | Algebraic Multiplicity | Geometric Multiplicity | Eigenspace Basis                          |
|------------|------------------------|------------------------|------------------------------------------|
| 4          | 2                      | 1                      | \(\{ \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} \}\) |
| 2          | 1                      | 1                      | \(\{ \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix} \}\) |

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
- Algebraic multiplicity is how many times an eigenvalue appears as a root of the characteristic polynomial.  
- Geometric multiplicity is the dimension of the eigenspace (number of independent eigenvectors) corresponding to that eigenvalue.  
- Here, eigenvalue 4 has algebraic multiplicity 2 but geometric multiplicity 1, meaning there is only one independent eigenvector for \(\lambda=4\).

---

### 5. Diagonalizability Conclusion

⚠️ Since the geometric multiplicity of \(\lambda = 4\) (which is 1) is **less than** its algebraic multiplicity (which is 2), the matrix \(A\) is **not diagonalizable**.

<span style="color:#800080; font-weight:bold;">💜 Final conclusion:</span> \(\boxed{A \text{ is not diagonalizable}}\).

<span style="color:#228B22; font-weight:bold;">💡 Explanation:</span>  
A matrix is diagonalizable if it has enough independent eigenvectors to form a basis — meaning the sum of geometric multiplicities equals the size of the matrix. Here, because the geometric multiplicity is less than the algebraic multiplicity for \(\lambda=4\), we don't have enough eigenvectors to diagonalize \(A\).

---

✅ **Summary:**

- The characteristic polynomial is \((4-\lambda)^2(2-\lambda)\).
- Eigenvalues: \(4\) (multiplicity 2), \(2\) (multiplicity 1).
- Eigenvectors for \(\lambda=4\) span a 1-dimensional space.
- Eigenvectors for \(\lambda=2\) span a 1-dimensional space.
- Because the eigenspace dimension for \(\lambda=4\) is less than its multiplicity, \(A\) is **not diagonalizable**.

---

<span style="color:#228B22; font-weight:bold;">💡 Tip:</span> For diagonalizability, the sum of geometric multiplicities must equal \(n\), the size of the matrix.

<!-- MathJax for LaTeX on GitHub Pages -->
<script>
  window.MathJax = {
    tex: {
      inlineMath: [['$','$'], ['\\(','\\)']],
      displayMath: [['$$','$$'], ['\\[','\\]']]
    },
    svg: { fontCache: 'global' }
  };
</script>
<script async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
