# Pipeline Validation & Testing

To test the nextflow pipeline integration, you can run a dry-run / stub run:

```bash
nextflow run main.nf -c tests/test.config -stub
```

Or run with the test config on real data:
```bash
nextflow run main.nf -c tests/test.config
```
