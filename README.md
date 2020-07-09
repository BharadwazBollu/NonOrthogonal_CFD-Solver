NonOrthogonal_FVM_SIMPLE_2D Whcih solves Non Dimensional Laminar Governing Equations.
Solver is robust enough to handle 70+ Non Orthogoanlity with good precision and uses volumetric interpolation for better accuracy.
Solver can handle even highly skewed grids.
While giving co-ordinates of your own mesh, make sure that increments in rows is increment in x-direction.
Solver contains preprocessing script to compute surface values and other parameters required for solver
Solver contains postprocessing script to view results and also optional script to export to .plt format which can be viewed in ParaView/TecPlot.

# Lid Driven Cavity

Grid = 41x41, Reynolds Number = 100

![lidDrivenMesh](https://user-images.githubusercontent.com/68074795/87076912-ff4a2980-c23f-11ea-830e-40782b578a3e.jpg)
![U_lidDrivenCenterline](https://user-images.githubusercontent.com/68074795/87076982-1ab53480-c240-11ea-909d-4c06be902219.jpg)
![V_lidDrivenCenterline](https://user-images.githubusercontent.com/68074795/87076991-1c7ef800-c240-11ea-8f52-3b1da35024fb.jpg)
![lidDriven_Ucontour](https://user-images.githubusercontent.com/68074795/87077059-31f42200-c240-11ea-85ed-d870f1a96aff.jpg)
![lidDriven_Vcontour](https://user-images.githubusercontent.com/68074795/87077066-33bde580-c240-11ea-9d3a-78cf1067addb.jpg)

# Pipe Bend Flow

Grid = 41x41, Reynolds Number = 100

![pipeBend_Mesh](https://user-images.githubusercontent.com/68074795/87077113-433d2e80-c240-11ea-9bf5-3f87f694b17d.jpg)
![pipeBend_Ucontour](https://user-images.githubusercontent.com/68074795/87077132-49330f80-c240-11ea-9fa2-728268c7a890.jpg)
![pipeBend_Vcontour](https://user-images.githubusercontent.com/68074795/87077138-49cba600-c240-11ea-9965-3343973f74cd.jpg)
![pipeBend_Quiver](https://user-images.githubusercontent.com/68074795/87077126-47694c00-c240-11ea-8a3f-b925a70782f4.jpg)


# Flow over Cylinder/Circle  

Grid = 41x41, Reynolds Number = 30

![Mesh](https://user-images.githubusercontent.com/68074795/87077210-6ec01900-c240-11ea-825b-cbd5afc73c0d.gif)
![bluffCircle_Ucontour](https://user-images.githubusercontent.com/68074795/87077227-74b5fa00-c240-11ea-9c96-1aa05da6ff4b.jpg)
![bluffCircle_Vcontour](https://user-images.githubusercontent.com/68074795/87077229-75e72700-c240-11ea-892c-9aae297aa501.jpg)
![bluffCircle_Quiver](https://user-images.githubusercontent.com/68074795/87077239-78498100-c240-11ea-879c-2751af173309.jpg)
