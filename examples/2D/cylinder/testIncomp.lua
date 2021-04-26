Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.1,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-0.1, -1.6, 5, 1.6},
		mshFile = "examples/2D/cylinder/geometry.msh"
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
		mu = 200,
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		WallsFixed = true,
		CylinderWallFixed = true,
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
			bodyForce = {0, 0},
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {1, 0, 0}
end

function Problem.IC:initWallsStates(pos)
	return {0, 0, 0}
end

function Problem.IC:initCylinderWallStates(pos)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:WallsV(pos, initPos, states, t)
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:CylinderWallV(pos, initPos, states, t)
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:FluidInputV(pos, initPos, states, t)
	return {1, 0}
end