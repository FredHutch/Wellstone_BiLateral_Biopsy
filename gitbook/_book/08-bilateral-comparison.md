# Bilateral comparison {#bilateral-comparison}





## STIR pairs simulation



## Fat infiltration
The Pearson correlation between right and left fat infiltration is 0.82.

```
#> # A tibble: 1 × 1
#>     cor
#>   <dbl>
#> 1 0.821
```

<img src="08-bilateral-comparison_files/figure-html/viz-fat-infilt-bilat-1.png" width="672" />

## Muscle strength
<img src="08-bilateral-comparison_files/figure-html/viz-muscle-strenth-bilat-1.png" width="672" />

## Comprehensive for other histopathological variables



<img src="08-bilateral-comparison_files/figure-html/viz-comprehensive-pathology-bilateral-1.png" width="672" />
## Baskets biolateral comparisons

### Per-gene correlation
We calculated the left-t-right correlation for every genes in the baskets. The unit of the expression here is TPM.
1. 


  
<table>
<caption>(\#tab:cor-and-mean-cor-per-gene-per-baksets)Mean of correlation of genes in each baskets calculated by TPM on L/R muscle.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> basket </th>
   <th style="text-align:right;"> cor_mean </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> DUX4-M6 </td>
   <td style="text-align:right;"> 0.6521460 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ECM </td>
   <td style="text-align:right;"> 0.5126943 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Inflamm </td>
   <td style="text-align:right;"> 0.4809635 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Complement </td>
   <td style="text-align:right;"> 0.6574381 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IG </td>
   <td style="text-align:right;"> 0.8281996 </td>
  </tr>
</tbody>
</table>

<img src="08-bilateral-comparison_files/figure-html/boxplot-per-gene-bilat-correlation-by-basket-1.png" width="672" />

__Scatter plot for each genes__
6 by 5 scatter plots
<img src="08-bilateral-comparison_files/figure-html/per-gene-correlation-scatter-DUX4-M6-1.png" width="672" />

<img src="08-bilateral-comparison_files/figure-html/per-gene-correlation-scatter-ECM-M6-1.png" width="672" />
<img src="08-bilateral-comparison_files/figure-html/per-gene-correlation-scatter-inflam-1.png" width="672" />
<img src="08-bilateral-comparison_files/figure-html/per-gene-correlation-scatter-Complement-1.png" width="672" />
<img src="08-bilateral-comparison_files/figure-html/per-gene-correlation-scatter-IG-1.png" width="672" />
### Per-basket by log sum
Here we calculate the per-basket bilateral correlation using basket scores.


<img src="08-bilateral-comparison_files/figure-html/viz-basket-score-bilateral-cor-1.png" width="672" />
