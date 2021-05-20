Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.006,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-0.01, -0.01, 0.628, 100},
		exclusionZones = {{0.29, 0, 0.316, 0.048}},
		mshFile = "examples/2D/damBreakWithObstacle/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "MinMax",
			outputFile = "tipPosition.txt",
			timeBetweenWriting = 0.01,
			minMax = "max",
			coordinate = 0 
		},
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.05,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		BoundaryFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.001,
		initialDT = 0.001,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, -9.81},
			computePres = false,
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end