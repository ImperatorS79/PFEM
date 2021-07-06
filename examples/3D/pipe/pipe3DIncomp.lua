Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 6,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.1,
		alpha = 1.2,
		omega = 0.37,
		gamma = 0.7,
		boundingBox = {-1, -1, -0.5, 5, 1, 1.5},
		addOnFS = true,
		deleteFlyingNodes = false,
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/3D/pipe/geometry.msh"
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
		maxDT = 0.001,
		initialDT = 0.001,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, 0, 0},
			gammaFS = 0.5,
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {1, 0, 0, 0}
end

function Problem.IC:initBoundaryStates(pos)
	return {0, 0, 0, 0}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, t)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:FluidInputV(pos, t)
	return {1, 0, 0}
end