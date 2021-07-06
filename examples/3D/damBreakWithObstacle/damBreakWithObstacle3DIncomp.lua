Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 1,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.1,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-1, -0.6, -0.01, 4, 0.6, 100},
		exclusionZones = {{2.396, -0.2, 0, 2.556, 0.2, 0.16}},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/3D/damBreakWithObstacle/geometry.msh"
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
		mu = 1e-3,
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
		maxDT = 0.01,
		initialDT = 0.01,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, 0, -9.81},
			gammaFS = 0.5,
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	local rho = Problem.Material.rho
	local g = -Problem.Solver.MomContEq.bodyForce[3]
	local z0 = 0.55
	
	if(pos[3] <= z0 and pos[1] <= 1.228 + 1.1*Problem.Mesh.hchar) then
		local p = rho*g*(z0 - pos[3])
		return {0, 0, 0, p}
	else 
		return {0, 0, 0, 0}
	end
end

function Problem.Solver.MomContEq.BC:WallsV(pos, t)
	return {0, 0, 0}
end