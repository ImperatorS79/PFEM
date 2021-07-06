Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 10.1,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.05,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		addOnFS = true,
		deleteFlyingNodes = false,
		exclusionZones = {},
		boundingBox = {-2, -2, -2, 2, 2, 2},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/3D/cubeToSphere/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.1,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		} 
	},
	
	Material = {
		mu = 1.0,
		rho = 100,
		gamma = 1.9
	},
	
	IC = {

	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.001,
		initialDT = 0.0000001,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, 0, 0},
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0, 0}
end
