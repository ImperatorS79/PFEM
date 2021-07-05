Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 1,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.025,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.55,
		boundingBox = {-0.05, -0.05, 10.05, 100},
		exclusionZones ={},
		addOnFS = true,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/ballInFluid/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.02,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		rho = 1000,
		mu = 1,
		gamma = 0.0
	},
	
	IC = {
		BoundaryFixed = true,
		DiskBoundaryFixed = false
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.01,
		initialDT = 0.01,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, -9.81},
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0}
end

function Problem.IC:initDiskBoundaryStates(pos)
	return {0, 1, 0}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, t)
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:DiskBoundaryV(pos, t)
	return {0, 1}
end