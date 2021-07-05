Problem = {
    id = "Bingham",
	simulationTime = 6,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.01,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		addOnFS = true,
		deleteFlyingNodes = true,
		boundingBox = {-0.01, -0.01, 0.99, 100},
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/dbDellaVecchia/geometry.msh"
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
			kind = "Mass",
			outputFile = "mass.txt",
			timeBetweenWriting = 0.01
		},
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.1,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 300,
		tau0 = 50,
		mReg = 1000,
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		WallsFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.005,
		initialDT = 0.000005,
		
		MomContEq = {
			minRes = 1e-4,
			maxIter = 10,
			gammaFS = 0,
			bodyForce = {0, -9.81},
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	local g = -Problem.Solver.MomContEq.bodyForce[2]
	local rho = Problem.Material.rho
	local p = 0
	if(pos[2] <= 0.267 and pos[1] <= 0.32 + 2*Problem.Mesh.hchar) then
		p = rho*g*(0.267 - pos[2])
	end
	return {0, 0, p}
end

function Problem.Solver.MomContEq.BC:WallsV(pos, t)
	return {0, 0}
end