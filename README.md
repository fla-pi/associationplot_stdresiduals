# Standardized residuals association plot
This script consists of five modified functions from <code>vcd</code> package (Meyer, Zeileis, Hornik 2021), and it allows to create association plots standardized residuals for chi-squared test.

It builds mainly on <code>assoc()</code> function, that returns residual-based shaded association plots from Pearson residuals.

## Modified Functions
The functions in the script are:

* <code>aspl()</code>: modified version of <code>assoc()</code>. This function returns the association plots.
* <code>struct()</code>: modified version of <code>strucplot()</code>.
* <code>legend_res()</code>: modified version of <code>legend_resbalance()</code>. This modified version was taken from [here](https://stackoverflow.com/questions/54674441/how-to-format-p-values-in-an-association-plot-made-with-vcd).
* <code>shading_hcl2()</code>: modified version of <code>shading_hcl()</code>.
* <code>vcdViewport()</code>: we kept this function as originally implemented in <code>vcd</code>.


## Usage
We use <code>Arthritis</code> dataset from <code>vcd</code> package, testing the association between <code>Treatment</code> and <code>Improved</code>.

```
> data("Arthritis")
> treat <- Arthritis$Treatment
> impr <- Arthritis$Improved
> table <- table(treat, impr)
> table
         impr
treat     None Some Marked
  Placebo   29    7      7
  Treated   13    7     21
  
> test2 <- chisq.test(treat, impr)

	Pearson's Chi-squared test

data:  treat and impr
X-squared = 13.055, df = 2, p-value = 0.001463

> test2$stdres
         impr
treat            None        Some      Marked
  Placebo  3.27419655 -0.09761768 -3.39563632
  Treated -3.27419655  0.09761768  3.39563632
```


First of all, you need to run the full script (you need to specify the path):

```
source("vcd_assoc_stdres.R", echo = TRUE)
```

Then, you can use the <code>aspl()</code> function. This function takes the contingency table as its only argument:

```
> aspl(table)
```

![Alt text](https://github.com/fla-pi/stdres_assocplot/blob/main/Rplot01.gif)


You can also rotate the graph by passing the optional argument <code>transform = T</code>:

```
> aspl(table, transform = T)
```

![Alt text](https://github.com/fla-pi/stdres_assocplot/blob/main/Rplot_transformed.gif)

Finally, it is also possible to create an association plot from the standardized residuals of Monte Carlo simulations. You can do this by passing the argument <code>simulate = T</code>. You can also change the number of iterations via the <code>B</code> argument (the default is 2000):

```
> aspl(table, simulate = T, B = 10000)
```

![Alt text](https://github.com/fla-pi/stdres_assocplot/blob/main/Rplot02.gif)

## Packages

Meyer, D., Zeileis, A., Hornik, K. (2021). _vcd: Visualizing Categorical Data_. R package version 1.4-9.
