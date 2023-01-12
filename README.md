T. Easley, Z. Ren, B. Kim, G. Karczmar, R.F. Barber, F.D. Pineda.

Nov. 2, 2020

# ECA Demo

ECA Demo is a MATLAB repository to share the code and a small-scale example of the simulator and reconstruction engine described in "Enhancement-Constrained Acceleration: A Robust Reconstruction Framework in Breast DCE-MRI" [1].

The minimal dataset of "Enhancement-Constrained Acceleration: A Robust Reconstruction Framework in Breast DCE-MRI" is publicly available and can be found [here](https://uchicago.app.box.com/s/3vfvdot8hzqausgr0mwfd85iwqce0dkq).

More detail on the reconstruction algorithm and input/output data strucutres can be found in **Algorithms_and_data_structures.pdf**.

## Installation

Zip and download ECA_demo and extract into a location on the MATLAB path. To keep only the simulation and reconstruction kernels, download the contents of **ECA_demo/Demo/demo_fctns**, ignoring **AP Phantom Functions** and **Paths**.

## Usage

```matlab
>> run_demo(1);  %% 0, 1, 2: demo sampling path options
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)

### References:

[1] T. O. Easley, Z. Ren, B. Kim, G. S. Karczmar, R. F. Barber, and F. D. Pineda, “Enhancement-constrained acceleration: A robust reconstruction framework in breast DCE-MRI,” PLOS ONE, vol. 16, no. 10, p. e0258621, Oct. 2021, doi: 10.1371/journal.pone.0258621.

 
### Associated Grants:

DMS-1654076 (National Science Foundation)

5R01CA218700-04 (National Institutes of Health)
