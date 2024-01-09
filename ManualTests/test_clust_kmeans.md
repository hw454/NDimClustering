# Description

The clustering program takes in GWAS data. This includes a matrix *β*
which contains the strength of association with a trait. The trait forms
the column and the snp is the row. There’s also a matrix of *p*-values
these represented the probability that the observed association is due
to chance. There’s also the *S**E* matrix with the standard error
association with the association scores for each snp, trait pair.

The snp’s are chosen based on those with the strongest association to
the exposure.

The algorithm attempts to use the association with different traits to
identify separate pathways within the data. We expect the pathways to
appear as lines in the original data-space, in order to identify
clusters we use two methods to reformat the data. The first is
transforming into the space given by the angles to the axes, this should
take lines through the origin with different gradients to balls of data.
The second is principal component analysis, which identifies to the
strongest trends within the data. In this testing process we will look
at different combinations of these.

For the test cases I have created synthetics which intentionally has
pathways in.

### Test Case Table

<table style="width:100%;">
<colgroup>
<col style="width: 14%" />
<col style="width: 14%" />
<col style="width: 14%" />
<col style="width: 14%" />
<col style="width: 14%" />
<col style="width: 14%" />
<col style="width: 14%" />
</colgroup>
<thead>
<tr class="header">
<th>Test No.</th>
<th>No. of pathways</th>
<th><a href="#pca">PCA method</a></th>
<th><a href="#angles">Angles used (Y or N)</a></th>
<th><a href="#cluster_p_weight">Weight cluster contribution with
p-values.</a></th>
<th>Clustering method</th>
<th><a href="#how_cents">Centroid assignment</a></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><a href="#test1">1</a></td>
<td>2</td>
<td>none</td>
<td>N</td>
<td>N</td>
<td>standard k-means</td>
<td>random</td>
</tr>
<tr class="even">
<td><a href="#test2">2</a></td>
<td><strong>4</strong></td>
<td>none</td>
<td>N</td>
<td>N</td>
<td>standard k-means</td>
<td>random</td>
</tr>
<tr class="odd">
<td><a href="#test3">3</a></td>
<td>2</td>
<td>none</td>
<td><strong>Y</strong></td>
<td>N</td>
<td>standard k-means</td>
<td>random</td>
</tr>
<tr class="even">
<td><a href="#test4">4</a></td>
<td><strong>4</strong></td>
<td>none</td>
<td>Y</td>
<td>N</td>
<td>standard k-means</td>
<td>random</td>
</tr>
<tr class="odd">
<td><a href="#test6">6</a></td>
<td>4</td>
<td>none</td>
<td>Y</td>
<td>N</td>
<td>standard k-means</td>
<td><strong>random from points</strong></td>
</tr>
<tr class="even">
<td><a href="#test7">7</a></td>
<td>4</td>
<td><strong>prcomp</strong></td>
<td>N</td>
<td>N</td>
<td>standard k-means</td>
<td>random from points</td>
</tr>
<tr class="odd">
<td><a href="#test8">8</a></td>
<td>4</td>
<td>prcomp</td>
<td><strong>Y</strong></td>
<td>N</td>
<td>standard k-means</td>
<td>random from points</td>
</tr>
<tr class="even">
<td><a href="#test9">9</a></td>
<td>4</td>
<td>prcomp</td>
<td>Y</td>
<td><strong>Y</strong></td>
<td>standard k-means</td>
<td>random from points</td>
</tr>
<tr class="odd">
<td><a href="#test10">10</a></td>
<td>4</td>
<td><strong>none</strong></td>
<td><strong>N</strong></td>
<td><strong>N</strong></td>
<td><a href="#min_aic"><strong>min-aic k-means</strong></a></td>
<td>random from points</td>
</tr>
<tr class="even">
<td><a href="#test11">11</a></td>
<td>4</td>
<td>none</td>
<td>Y</td>
<td>N</td>
<td><strong>min-aic k-means</strong></td>
<td>random from points</td>
</tr>
<tr class="odd">
<td><a href="#test12">12</a></td>
<td>4</td>
<td>none</td>
<td>Y</td>
<td><strong>Y</strong></td>
<td>min-aic k-means</td>
<td>random from points</td>
</tr>
<tr class="even">
<td><a href="#test13">13</a></td>
<td>4</td>
<td><strong>prcomp</strong></td>
<td><strong>N</strong></td>
<td><strong>N</strong></td>
<td>min-aic k-means</td>
<td>random from points</td>
</tr>
<tr class="odd">
<td><a href="#test14">14</a></td>
<td>4</td>
<td>prcomp</td>
<td><strong>Y</strong></td>
<td>N</td>
<td>min-aic k-means</td>
<td>random from points</td>
</tr>
<tr class="even">
<td><a href="#test15">15</a></td>
<td>4</td>
<td>prcomp</td>
<td>Y</td>
<td><strong>Y</strong></td>
<td>min-aic k-means</td>
<td>random from points</td>
</tr>
<tr class="odd">
<td><a href="#test16">10</a></td>
<td>4</td>
<td><strong>none</strong></td>
<td><strong>N</strong></td>
<td><strong>N</strong></td>
<td><a href="#dbscan"><strong>dbscan</strong></a></td>
<td>centroids not assigned in density based methods</td>
</tr>
<tr class="even">
<td><a href="#test17">11</a></td>
<td>4</td>
<td>none</td>
<td>Y</td>
<td>N</td>
<td><strong>dbscan</strong></td>
<td>centroids not assigned in density based methods</td>
</tr>
<tr class="odd">
<td><a href="#test18">13</a></td>
<td>4</td>
<td><strong>prcomp</strong></td>
<td><strong>N</strong></td>
<td><strong>N</strong></td>
<td>dbscan</td>
<td>centroids not assigned in density based methods</td>
</tr>
<tr class="even">
<td><a href="#test19">14</a></td>
<td>4</td>
<td>prcomp</td>
<td><strong>Y</strong></td>
<td>N</td>
<td>dbscan</td>
<td>centroids not assigned in density based methods</td>
</tr>
<tr class="odd">
<td><a href="#test20">15</a></td>
<td>4</td>
<td>prcomp</td>
<td>Y</td>
<td>N</td>
<td>dbscan</td>
<td>centroids not assigned in density based methods</td>
</tr>
<tr class="even">
<td><a href="#test20">15</a></td>
<td>4</td>
<td>prcomp</td>
<td>Y</td>
<td>N</td>
<td>dbscan</td>
<td>centroids not assigned in density based methods</td>
</tr>
</tbody>
</table>

# Instructions for using the package.

## Install package

Load the “ndimclusteringR” package

## Iteration parameters

The parameters for the program are:

-   `dname` The directory where the data is saved. This should include
    the following files:

    -   `unstdBeta_df.csv` A `.csv` file with columns seperated by `,`
        and rows by `;` , each row corresponds to a snp and each column
        is a trait. The values correspond to the intensity of
        association determined from a GWAS. The data is sourced from
        openGWAS.

    -   `unstdSE_df.csv` A `.csv` file with columns seperated by `,` and
        rows by `;` , each row corresponds to a snp and each column is a
        trait. The values correspond to the standard error of the
        association determined from a GWAS. The data is sourced from
        openGWAS.

    -   `pval_df.csv` A `.csv` file with columns seperated by `,` and
        rows by `;` , each row corresponds to a snp and each column is a
        trait. The values correspond to the *p*-values of the
        association determined from a GWAS. The *p*-values are the
        probability that the intensity values are observed due to
        chance. The data is sourced from openGWAS.

    -   `trait_info_nfil.csv` A `.csv` file with columns seperated by
        `,` and rows by `;` , each row corresponds to a trait. The label
        for the trait is `phenotype`, to get the full definition for the
        trait as used in the data collection this is found in the column
        `description`.

-   `clust_type` method of clustering to be used.

-   `basic` is a standard *k*-means clustering method.

-   `min` is an iteration of *k*-means, which minimises the *a**i**c*
    from cluster sets run from 1 cluster to *k*.

-   `dbscan` density based clustering method.

-   `pca_type` method for determining the principle components. See [PCA
    section](#pca)

    -   `prcomp` PCA is performed using the `prcomp` package.
    -   `none` no PCA is performed and the data is clustered without
        PCA.

-   `nclust` The number of clusters to search for when using `basic` or
    `min` clustering.

-   `n_pcs` The number of principal components to reduce to. This term
    is only used if `pca_type` is NOT `none`.

-   `bin_angles` Binary switch determining whether the data should be
    transformed to the angle space or not. If 0 then nothing is done, if
    1 then the data is transformed from the intensity scores to the
    angle of those scores to the axis in that data-space. Since the
    angle describe a position betwen two axes by using the angles we can
    reduce the dimension of the data-space by 1, if we wanted to
    retained all information we would also store the distance from the
    origin but we are searching for elements on the same line, it does
    not matter to our clusters where on the line they are, therefore we
    do not store this distance.

-   `bin_p_clust` Binary switch. If 0 then *k*-means clustering
    calculated the cluster centroids based on the average of the terms
    in that cluster. If 1 then the cluster centroids are recalculated
    based on a weighted average using the intensity and the *p*-values.
    See [*p*-value weighting section](#cluster_p_weight).

-   `bin_p_score` Binary switch indicating the weighting when scoring
    clusters on trait axis. If 0 then the clusters are scored based on
    their elements average intensity scores on each axis. If 1 then the
    clusters are scored based on a weight average of their elements
    scores on each trait, weighted by the *p*-value that the intensity
    score is a true representation.

-   `how_cents` The method to be used for initialising the centres of
    the clusters.

    -   `rand` The value on each axis is randomly assigned using a
        uniform distribution between the max and min on each axis.
    -   `points` The centre for each cluster is assigned to points from
        the dataset. The points are randomly selected with a uniform
        distribution.

### What happens within the program? - Algorithm: Step by step

If all the steps are included then they are performed in the following
order.

1.  Crop data to complete cases.
2.  Convert data to angles.
3.  Perform PCA on angle data.
4.  Cluster on the PCA

-   Set cluster centres
-   Converge cluster centres
-   if cluster type “min” then rerun for different number of clusters
    and choose set with min-aic.

### Test 1

Run the tests with:

-   No PCA
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 2.

<!-- -->

    ## [1] "Clusters converged 3"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-3-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-3-2.png)

### Test 2

For test2 we will keep the parameters the same as [test1](#test1) but
increase the number of pathways.

Run the tests with:

-   No PCA
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   no angles calculated before clustering.
-   **number of pathways 4.**
-   Clusters centres are not weighted by the p-values.

<!-- -->

    ## [1] "Clusters converged 5"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-5-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-5-2.png)

## Calculating Angles

The data-points for each snp have an intensity score associated with
each trait. We want to identify lines within our data but we can see
from the first two tests that standard clustering does not detect these
lines. To try and convert data that is visualised as lines to something
that looks like balls of data we convert the data to the angle to each
axis.

1.  Create a unit-vector for each axis (you can remove one axis since
    the angles reduce the dimension), the algorithm removes the last
    axis.
2.  For each unit vector found and each data-point *x*:
    -   
        $$\theta = \frac{x \cdot u}{||u|| ||x||}$$
    -   If *θ* &gt; *π* reassign the angle *θ* = *θ* − *π*. These are
        angles on opposite sides of the circle but they lie on the same
        line so we want the algorithm to group them together.

### Test 3

For test3 we will return to the case of two pathways but this time we
will reformat the data before clustering. We will transform the data
from the raw *β* intensity scores to the angles of the scores when
transformed to a polar-coordinate system.

Run the tests with:

-   No PCA
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   **angles calculated.**
-   number of pathways 2.

<!-- -->

    ## [1] "Clusters converged 7"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-7-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-7-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-7-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-7-4.png)

### Test 4

For test4 we use the angles again but we investigate the case with 4
pathways.

Run the tests with:

-   No PCA
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   angles calculated before clustering.
-   **number of pathways 4.**

<!-- -->

    ## [1] "Clusters converged 4"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-9-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-9-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-9-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-9-4.png)

## Consider the location of the centroids.

The previous tests used centroids assigned randomly within the ranges on
each axis. The results show that although we ask for *k* clusters we
often get some of the clusters returning as empty. This is likely to be
caused by the initial centroids being outside the data-space.

Instead we will initialise the centroids using randomly selected points
from out dataset.

### Test 5

For test5 we keep the parameters the same as [test2](#test2) but we
change the initialisation of the cluster centroids.

Run the tests with:

-   No PCA
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.
-   **Centroids assigned using points.**

<!-- -->

    ## [1] "Clusters converged 9"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-11-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-11-2.png)

### Test 6

For test6 we keep the parameters the same as [test5](#test5) but we
include the conversion to angles.

Run the tests with:

-   No PCA
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   **angles calculated before clustering.**
-   number of pathways 4.
-   Centroids assigned using points

<!-- -->

    ## [1] "Clusters converged 3"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-13-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-13-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-13-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-13-4.png)

## Principal Components

So far we have ignored the principal components option. This is another
method for reducing the dimension. It identifies the eigen-vectors and
corresponding eigen-values of the data space. Then outputs the set of
eigen-vectors with the largest eigen-values as a matrix. This matrix is
then used to transform the data onto a new data-space with reduced
dimension. The principal components characterise the variation in the
data. By defining the axes as the eigen-vectors we remove the highest
levels of correlation within the data. In doing this we are trying to
find the level at which the data separates into distinct categories.

PCA steps:

1.  Ensure each data axis is normalised.
2.  Identify the eigen-vectors and eigen-values of the data-space.
3.  Use the eigen-vectors of the `n_pc` largest eigen-values to form the
    columns of a matrix.
4.  Take each data-point as a row-vector and multiply by the matrix to
    transform the data onto a new data-space.
5.  Transform the p-value matrix and standard error matrix onto the same
    dataspace.
6.  Store the transform matrix. This represents the proportion of each
    trait which makes up each pc vector.

[`prcomp`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp)
is a pca package in R.

### Test 7

For the first pca test we will turn off the calculation of angles and
use a standard k-means clustering to try and visualise what the pca is
doing on it’s own. This is [test5](#test5) with `how_cents` changed.

Run the tests with:

-   **PCA method prcomp**
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points

<!-- -->

    ## [1] "Clusters converged 5"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-15-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-15-2.png)

### Test 8

For the second pca test we will include the calculation of angles and
use a standard k-means clustering to try and visualise what how the pca
works with the angle data. This is [test6](#test6) with `how_cents`
changed.

Run the tests with:

-   PCA method prcomp
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   **angles calculated before clustering.**
-   number of pathways 4.
-   Centroids assigned using points

<!-- -->

    ## [1] "Clusters converged 3"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-17-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-17-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-17-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-17-4.png)

## Weight cluster centroid by p-value

Each data-point has an associated p-value. This corresponds to the
chance the score observed is due to chance. We do not want to skew our
conclusions by unreliable snps. Therefore when the cluster centroids are
found we will use a weighted average based on the p-value. For each
co-ordinate in the cluster centroid the new centroid after cluster
assignment is given by:

$$
\pmb{c}\_i = \frac{\sum\_{\text{points in cluster}}\beta\_i(1-p\_i)}{\sum\_{\text{points in cluster}}(1-p\_i)}
$$

### Test 9

Run the tests with:

-   PCA method prcomp
-   basic k-means clustering.
    -   k is set to the number of pathways.
-   angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points
-   **Centroids reassigned during clustering using a p-value
    weighting.**

<!-- -->

    ## [1] "Clusters converged 1"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-19-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-19-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-19-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-19-4.png)

## Min-AIC

So far we have clustered using a fixed number of clusters. Since we are
using synthetic test data we can set the value for *k* to the number of
clusters we expect. However when data-mining we don’t know how many
clusters we are looking for. We consider an iterative clustering
approach, running a standard k-means clustering on 1 to *k* clusters.
For each set of clusters calculate the Akaike Information Criterion
(AIC), then proceed with the cluster set that minimises this.

The AIC is given by:
$$ 
l = \frac{\sum\_{\text{data\_points}}\text{dist}\left(\beta\_i, c\_{\text{centre for cluster assigned to } \beta\_i}\right)^2}{\text{var}(\text{dist}\left(\pmb{\beta}\right))}, \\
d = \text{dimension of each point},\\
k = \text{number of points}, \\
AIC = d + 2  k  l
$$

For *n**c* in \[1,*k*\]

1.  Set the number of clusters to *n**c*.
    -   If *n**c* &gt; *k* go to step 6.
2.  Run k-means clustering for *k* = *n**c*.
3.  Calculate the *a**i**c* for the clusters.
4.  Increase *n**c* by 1.
5.  Go to step 1.
6.  Find the number of clusters which gives the minimum aic.
7.  Return the cluster set for this number of clusters.

We will perform the following tests on the min-aic clustering method:

<table>
<colgroup>
<col style="width: 13%" />
<col style="width: 13%" />
<col style="width: 13%" />
<col style="width: 13%" />
<col style="width: 17%" />
<col style="width: 13%" />
<col style="width: 13%" />
</colgroup>
<thead>
<tr class="header">
<th>Test No.</th>
<th>No. of pathways</th>
<th><a href="#pca">PCA method</a></th>
<th><a href="#angles">Angles used (Y or N)</a></th>
<th><a href="#cluster_p_weight">Weight cluster contribution with
p-values.</a></th>
<th>Clustering method</th>
<th><a href="#how_cents">Centroid assignment</a></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><a href="#test10">10</a></td>
<td>4</td>
<td><strong>none</strong></td>
<td><strong>N</strong></td>
<td><strong>N</strong></td>
<td><a href="#min_aic"><strong>min-aic k-means</strong></a></td>
<td>random from points</td>
</tr>
<tr class="even">
<td><a href="#test11">11</a></td>
<td>4</td>
<td>none</td>
<td>Y</td>
<td>N</td>
<td><strong>min-aic k-means</strong></td>
<td>random from points</td>
</tr>
<tr class="odd">
<td><a href="#test12">12</a></td>
<td>4</td>
<td>none</td>
<td>Y</td>
<td><strong>Y</strong></td>
<td>min-aic k-means</td>
<td>random from points</td>
</tr>
<tr class="even">
<td><a href="#test13">13</a></td>
<td>4</td>
<td><strong>prcomp</strong></td>
<td><strong>N</strong></td>
<td><strong>N</strong></td>
<td>min-aic k-means</td>
<td>random from points</td>
</tr>
<tr class="odd">
<td><a href="#test14">14</a></td>
<td>4</td>
<td>prcomp</td>
<td><strong>Y</strong></td>
<td>N</td>
<td>min-aic k-means</td>
<td>random from points</td>
</tr>
<tr class="even">
<td><a href="#test15">15</a></td>
<td>4</td>
<td>prcomp</td>
<td>Y</td>
<td><strong>Y</strong></td>
<td>min-aic k-means</td>
<td>random from points</td>
</tr>
</tbody>
</table>

### Test 10

Run the tests with:

-   No PCA
-   **min-aic k-means clustering.**
    -   Max *k* is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points
-   Centroids are not reassigned during clustering using a p-value
    weighting.

<!-- -->

    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 3"
    ## [1] "Clusters converged 6"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-21-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-21-2.png)

### Test 11

For test 11 we keep everything the same as [test10](#test10) with the
addition of the conversion to angles before clustering.

Run the tests with:

-   No PCA
-   min-aic k-means clustering.
    -   Max *k* is set to the number of pathways.
-   **angles calculated before clustering.**
-   number of pathways 4.
-   Centroids assigned using points
-   Centroids are not reassigned during clustering using a p-value
    weighting.

<!-- -->

    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 3"
    ## [1] "Clusters converged 3"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-23-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-23-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-23-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-23-4.png)

### Test 12

For test 12 we keep everything the same as [test11](#test11) with the
addition of weighting the cluster centre’s by the p-values.

Run the tests with:

-   No PCA
-   min-aic k-means clustering.
    -   Max *k* is set to the number of pathways.
-   angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points
-   **Centroids are reassigned during clustering using a p-value
    weighting.**

<!-- -->

    ## [1] "Clusters converged 1"
    ## [1] "Clusters converged 1"
    ## [1] "Clusters converged 1"
    ## [1] "Clusters converged 1"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-25-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-25-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-25-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-25-4.png)

### Test 13

For test 13 we return to the set up of the first min-aic test in
[test10](#test10) with the addition of principal component analysis.

Run the tests with:

-   **PCA using prcomp**
-   min-aic k-means clustering.
    -   Max *k* is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points
-   Centroids are not reassigned during clustering using a p-value
    weighting.

<!-- -->

    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 4"
    ## [1] "Clusters converged 5"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-27-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-27-2.png)

### Test 14

For test 14 we add the angle conversion to [test13](#test13).

Run the tests with:

-   PCA using prcomp
-   min-aic k-means clustering.
    -   Max *k* is set to the number of pathways.
-   **angles calculated before clustering.**
-   number of pathways 4.
-   Centroids assigned using points
-   Centroids are not reassigned during clustering using a p-value
    weighting.

<!-- -->

    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 2"
    ## [1] "Clusters converged 4"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-29-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-29-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-29-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-29-4.png)

### Test 15

For test 15 we take the set up in [test14](#test14) and we [weight the
cluster centre’s using the p-value weighting](#cluster_p_weight).

Run the tests with:

-   PCA using prcomp
-   min-aic k-means clustering.
    -   Max *k* is set to the number of pathways.
-   angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points
-   [**Centroids are reassigned during clustering using a p-value
    weighting.**](#cluster_p_weight)

<!-- -->

    ## [1] "Clusters converged 1"
    ## [1] "Clusters converged 1"
    ## [1] "Clusters converged 1"
    ## [1] "Clusters converged 1"

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-31-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-31-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-31-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-31-4.png)

## DBSCAN

**dbscan** is a density based clustering method. Points are assigned to
the cluster if there are within range of another point in the cluster.
As a result the clusters can be different shapes and not necessarily
balls.

Dbscan does not require a predefined number of clusters, it instead
takes in a set density for the points `point_eps`.

We will repeat the tests we used for [min-aic](#minaic) changing the
clustering method to dbscan, however we won’t vary the centroid
calculations as centroids aren’t used in density clustering.

### Test 16

Run the tests with:

-   No PCA
-   **dbscan clustering.**
-   no angles calculated before clustering.
-   number of pathways 4.
-   Centroids not assigned in density based methods.

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-33-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-33-2.png)

### Test 17

Take the case from [test 16](#test16) and introduce the conversion to
angles before clustering.

Run the tests with:

-   No PCA
-   dbscan clustering.
-   **Angles calculated before clustering.**
-   number of pathways 4.
-   Centroids not assigned in density based methods.

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-35-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-35-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-35-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-35-4.png)

### Test 18

Take the case from [test 16](#test16) and introduce the conversion to
principal components before clustering. First we will do this without
calculating the angles.

Run the tests with:

-   No PCA
-   dbscan clustering.
-   **Angles calculated before clustering.**
-   number of pathways 4.
-   Centroids not assigned in density based methods.

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-37-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-37-2.png)

### Test 19

Now let’s combine [test 17](#test17) and [test 18](#test18) and
calculate the principal components and the angles before clustering.

Run the tests with:

-   **PCA calculated using prcomp**
-   dbscan clustering.
-   **Angles calculated before clustering.**
-   number of pathways 4.
-   Centroids not assigned in density based methods.

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-39-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-39-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-39-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-39-4.png)

## Vary the density

### Test 20

Now let’s investigate the density parameter to **dbscan**. Let’s make it
smaller and try `point_eps = 0.2`.

Run the tests with:

-   No PCA
-   dbscan clustering.
-   **Angles calculated before clustering.**
-   number of pathways 4.
-   Centroids not assigned in density based methods.

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-41-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-41-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-41-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-41-4.png)

### Test 21

Let’s test the density parameter again but this time increase to
`point_eps = 0.8`

Run the tests with:

-   No PCA
-   dbscan clustering.
-   **Angles calculated before clustering.**
-   number of pathways 4.
-   Centroids not assigned in density based methods.

![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-43-1.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-43-2.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-43-3.png)![](test_clust_kmeans_files/figure-markdown_strict/unnamed-chunk-43-4.png)
