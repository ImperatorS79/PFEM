Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.05,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.3,
		boundingBox = {-2, -2, 2, 2},
		exclusionZones = {},
		addOnFS = true,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/squareToDisk/geometry.msh"
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
		maxDT = 0.0025,
		initialDT = 0.0025,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, 0},
			gammaFS = 1,
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0}
end
