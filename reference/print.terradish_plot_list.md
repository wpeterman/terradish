# Print terradish panel plot lists

Internal helper used so multi-panel terradish plot objects auto-print to
the active graphics device without printing list structure to the
console.

## Usage

``` r
# S3 method for class 'terradish_plot_list'
print(x, ...)
```

## Arguments

- x:

  A `terradish_plot_list` object.

- ...:

  Unused.

## Value

The input `x`, a `terradish_plot_list` containing named `ggplot` panels,
returned invisibly. Each panel is printed to the active graphics device
for its graphical side effect.
