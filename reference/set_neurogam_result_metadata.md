# Attach analysis metadata without altering result fields or S3 classes

Stores versioned metadata as the \`neurogam_metadata\` attribute. This
keeps legacy list field names and classes unchanged and survives
\`saveRDS()\`. Existing metadata are replaced by fields supplied in the
new metadata list.

## Usage

``` r
set_neurogam_result_metadata(x, metadata)
```

## Arguments

- x:

  A recognized neurogam result.

- metadata:

  Named metadata returned by \`make_neurogam_metadata()\`.

## Value

The same result with a metadata attribute.
