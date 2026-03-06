# Developer notes
## Non-unique species names

According to the NCBI, novel unconfirmed species are no longer given unique names but a gathered under a single name bin,
e.g. `Bacillus sp.`.
It's not clear if they also share the ID.

```
Previously: Bacillus sp. St12345 
Now: Bacillus sp.
This means that in these cases sequences of several species will be mapping to a single node in the taxonomy database.
```

