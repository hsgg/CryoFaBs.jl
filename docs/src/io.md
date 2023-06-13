# IO

Calculating a *cryofab* can be time-consuming. Therefore, it may be desirable
to save it to disk. This can be done with the `write()` function, which works
with both 2D `AngularCryoFaB`s and 3D `AngRadCryoFaB`s.

Suppose `cfb` is an `CryoFaB`. Then, you can write to the file `filename.cfb` with
```julia
write("filename.cfb", cfb)
```
Here we used the `.cfb` file extension. Obviously, anything your operating
system supports works just as well, like `.bin`.

To read, we need to distinguish between an `AngularCryoFaB` and an
`AngRadCryoFaB`, e.g.,
```julia
cfb = AngularCryoFaB("filename.cfb")

cfb = AngRadCryoFaB("filename.cfb")
```


### Another way to read

There is also a more basic-looking way to read a *cryofab*:
```julia
cfb = read("filename.cfb", AngularCryoFaB)
```
for 2D, or
```julia
cfb = read("filename.cfb", AngRadCryoFaB)
```
for 3D.
