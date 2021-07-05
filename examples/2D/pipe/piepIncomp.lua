Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 6,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.05,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		addOnFS = true,
		deleteFlyingNodes = false,
		exclusionZones = {},
		boundingBox = {-0.01, -0.01, 5, 1.01},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/pipe/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.5,
			whatToWrite = {"p", "u"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 200,
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		BoundaryFixed = false,
		FluidInputFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.025,
		initialDT = 0.025,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, 0},
			gammaFS = 0.5,
			computePres = false,
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0}
end

function Problem.IC:initFluidInputStates(pos)
	return {1, 0, 0}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, t)
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:FluidInputV(pos, t)
	return {1, 0}
end