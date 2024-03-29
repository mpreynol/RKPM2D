# RKPM2D: 2D implementation of RKPM written by Mathew Reynolds in December of 2016. 

Code Structure:
   Included Classes:
      -QuadMesh: Defines background mesh. Uses a 2D FEM mesh with shape functions. This provides a robust mechanism to retrieve the spatial cordinates and jacobian value for an arbitrary isoparametric element.
      -GridRectangle: Defines the Nodes and Connectivity for a 2D rectangular mesh.
      -Boundary: Assigns nodes as a boundary value in an output array.
      -Quadrature: Quadrature tables in 1D and 2D. 2D tables are organized as the tensor product to reduce loops later.
      -Cloud: Defines a cloud of RKPM Nodes.
      -RKNode: Defines a class for each RKNodes, stores cordinates and values of shape and weight functions
      -RKShape: Defines the shape function for each RKNode.
      -Weight: Defines the kernal used for each node - note that the kernal for any paticular node may be singular.
      -Cloud2: Defines a cloud of RKPM Nodes, the domain and boundary integration subroutines are written to incorporate multithreading.
      -MeshPlot: Misc functions used to plot the background mesh. Most of the functionality is saved from FEM implementation.
      -Element: Class defined for each element of quadmesh, has functions to return spatial cordinates of guass quadrature points and jacobians.
      -Domains: A class that holds collections of "QuadMesh" for domain decomposition. Note that the QuadMesh runtimes save local copies of RK values.

Program Workflow (using test.m):
1. Background Mesh is defined as instance of "QuadMesh" or "Domains". The order of integration is specified globally, as well as the number of threads.
2. A RKPM point cloud is created as instance of "Cloud". This is done by using the 'Nodes' array from a Rectangular Grid or a random method writtin in the cloud objects.
3. A Check of Partial Unity, Partial Nullity, and Linear Consistency can be performed on an instance of QuadMesh.
4. Weight and RKPM shape function values are precomputed and stored in Element Instances. If the Mesh is:
	- A consistent instance of QuadMesh: Then the Stored Values can be exported directly as arrays to be sent by value into the integration routine.
	- An instance of "Domains": Then the stored values of each "QuadMesh" within in "Domains" is stored for when we call parallel integration.
5. The consituitive array is processed.
6. The Domain in integrated by calling the appropriate method within the instances of PointCloud. This either loops through elements, then quadrature points, then nodes*2; or, loops through the domains, then elements, then quadrature points, then nodes*2.  
7. The boundary is integrated either exactly if a function handle is passed into the arrays, or by linearly interpolating between assigned values of nodes in the BN array.
8. The system is then solved using a direct solve algorithm in matlab.

Remarks:
1. BE array defines the essential boundary array, this is used passed in "RKNodes" to determine if the node is singular or not (if !=-Inf).
2. BN array defines the natural boundary array, this is passed into the integration scheme to either: flag the correct edge for exact integration, or to perform linear interpolation on the actual values.
3. Domain Decomposition can be "turned off" by commenting and uncommeting several lines.
4. Domain Decomposition assigns bins of elements according to their order in the "NEL" array. Because of this the Decomposition may include non-contiguous elements.
5. Each Domain in the decomposition includes the possibility of having every node within it, this in the parallel loop that constructs the stiffness and force vector for each domain the full system array is worked on. This is a consequence of: Remark4, and I haven't written a better decomposition algorithm yet. As it stands this one should double the speed of integration using 4 cores vs. 1.
6. The arrays from each domain integration is summed in the '3rd' direction to form the complete array.
7. 'return Interpolated' does not use precomputed values of derivatives because the intent is that user can use a point that need not be a subset of the quadrauture points. This method is also used to plot the mesh deformation (currently weights and shapefunction values are not stored at the RKNode positions.
8. The random node generator that exist currently can only populate points from the origin out to mesh limits. You may add specified points (for imposiing essential bounday conditions) by passing in an array in the 'createPoints' module


