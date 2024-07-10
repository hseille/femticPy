# femticPy



`femticpy` is a Python module to create inputs / read outputs of  FEMTIC, a 3-D magnetotelluric inversion code (Usui, 2015).

`femticpy` is a python module that allows to: 

- create the input observed data file from a list of .edi files
- prepare the files required to generate a mesh
- read, plot and QC the result of the inversion 

## CSIRO Collection
[http://hdl.handle.net/102.100.100/636092?index=1](http://hdl.handle.net/102.100.100/636092?index=1)


## Workflow


The workflow consists in:
1. Create topography, bathymetry and coast line files using the notebook `prep_topo_bathy.ipynb`
2. Create input files to create the mesh and inversion files using the notebook `prep_inversion.ipynb`
3. Create the mesh using the pipeline `create_mesh.sh`, that calls the Pre/post-processing tools (available on FEMTIC Repo)
4. Run the inversion using `run_inversion.sh`, that calls the `femtic` executable
5. Analyse and plot the inversion results using the notebook `femtic_results.ipynb`


The module can be used using the Jupyter notebooks provided in the examples directory. 3 examples are provided: 
- A real data example in examples/auslamp, to create a mesh in a coastal area, run the inversion and QC inversion results (.edi files can be found on the Geological Survey of South Australia website:  [https://www.energymining.sa.gov.au/industry/geological-survey/gssa-projects/em-geophysics](https://www.energymining.sa.gov.au/industry/geological-survey/gssa-projects/em-geophysics))
- a synthetic MMT data example in examples/mmt, to create a mesh in a full marine environment using real bathymetry data
- a synthetic land MT example examples/synthetic, toy example that can be run quickly to test the inversion of Z and VTF data



## Notes

- Models can be visualized in Paraview
- The directory structure for each new project has to be copied from the `empty_project` directory
- The edi files read well when exported from winglink or geotools, other formats can cause problems
- The following codes must be compiled:
	- femtic (inversion executable)
	- makeMtr (Pre/post-processing tools)
	- makeTetraMesh (Pre/post-processing tools)
	- tetgen ([https://github.com/libigl/tetgen](https://github.com/libigl/tetgen))


## Dependencies

- pyproj
- rdp
- shapely
- pandas > 1.5


## References: 

### Papers:

1. Y. Usui, 3-D inversion of magnetotelluric data using unstructured tetrahedral elements: applicability to data affected by topography, Geophys. J. Int., 202 (2): 828-849, [https://doi.org/10.1093/gji/ggv186](https://doi.org/10.1093/gji/ggv186), 2015).

2. Y. Usui, Y. Ogawa, K. Aizawa, W. Kanda, T. Hashimoto, T. Koyama, Y. Yamaya and T. Kagiyama, Three-dimensional resistivity structure of Asama Volcano revealed by data-space magnetotelluric inversion using unstructured tetrahedral elements, Geophys. J. Int., 208 (3): 1359-1372, [https://doi.org/10.1093/gji/ggw459](https://doi.org/10.1093/gji/ggw459), 2017.

3. Y. Usui, T. Kasaya, Y. Ogawa and H. Iwamoto, Marine magnetotelluric inversion with an unstructured tetrahedral mesh, Geophys. J. Int., 214(2): 952-974, [https://doi.org/10.1093/gji/ggy171](https://doi.org/10.1093/gji/ggy171), 2018.


### FEMTIC repository:

[https://github.com/yoshiya-usui/femtic](https://github.com/yoshiya-usui/femtic)


### FEMTIC website:

[https://sites.google.com/view/yoshiyausui/femtic](https://sites.google.com/view/yoshiyausui/femtic)


### Resources / Manuals:

Manual_Of_FEMTIC_v4.1.pdf ([https://github.com/yoshiya-usui/femtic/tree/main/doc](https://github.com/yoshiya-usui/femtic/tree/main/doc))

How_To_Make_Tetra_Mesh_For_FEMTIC.pdf ([https://github.com/yoshiya-usui/makeTetraMesh/tree/master/doc](https://github.com/yoshiya-usui/makeTetraMesh/tree/master/doc))





## License

[MIT](https://choosealicense.com/licenses/mit/)