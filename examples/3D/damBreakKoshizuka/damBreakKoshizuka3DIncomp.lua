Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.0146,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.5,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-1, -1, -1, 1.584, 1.175, 100},
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/3D/damBreakKoshizuka/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.2,
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
		maxDT = 0.00025,
		initialDT = 0.00025,
		
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
	local g = -Problem.Solver.MomContEq.bodyForce[3]
	local rho = Problem.Material.rho
	local p = 0
	if(pos[3] <= 2*0.146 and pos[1] <= 0.146 + 2*Problem.Mesh.hchar) then
		p = rho*g*(2*0.146 - pos[3])
	end
	return {0, 0, 0, p}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, t)
	return {0, 0, 0}
end