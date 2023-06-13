# IO

Calculating a *cryofab* can be time-consuming. Therefore, it may be desirable
to save it to disk. This can be done with the `write()` function. Suppose `cfb`
is an `AngularCryoFaB`. Then, you can read and write as
```julia
write("filename.cfb", cfb)

cfb = AngularCryoFaB("filename.cfb")
```
`filename.cfb` is the filename it gets written to.

If it is a 3D `AngRadCryoFaB`, then, similarly,
```julia
write("filename.cfb", cfb)

cfb = AngRadCryoFaB("filename.cfb")
```

There is also a more basic-looking way to read a *cryofab*:
```julia
cfb = read("filename.cfb", AngularCryoFaB)
```
for 2D, or
```julia
cfb = read("filename.cfb", AngRadCryoFaB)
```
for 3D.
