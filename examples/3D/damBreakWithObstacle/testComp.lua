Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 1,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.02,
		alpha = 1.25,
		omega = 0.6,
		gamma = 0.35,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-0.1, -0.6, -0.1, 4, 0.6, 100},
		exclusionZones = {{2.396, -0.2, 0, 2.556, 0.2, 0.16}},
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
		gamma = 0,
		K0 = 220000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		WallsFixed = true,
	},
	
	Solver = {
	    id = "CDS_Meduri",
		adaptDT = true,
		maxDT = 0.001,
		maxRemeshDT = 0.05,
		initialDT = 1e-8,
		securityCoeff = 100,
		
		MomEq = {
			bodyForce = {0, 0, -9.81},
			BC = {
			
			}
		},
		
		ContEq = {
			strongContinuity = true,
			enableStab = true,
			bodyForce = {0, 0, -9.81},
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
    local rhoStar = Problem.Material.rhoStar
    local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local g = -Problem.Solver.MomEq.bodyForce[3]
	local z0 = 0.55
	
	if(pos[3] <= z0 and pos[1] <= 1.228 + 1.1*Problem.Mesh.hchar) then
	    local rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[3]) + 1)^(1/(K0p - 1))
		local p = K0/K0p*((rho/rhoStar)^K0p - 1)
		return {0, 0, 0, p, rho, 0, 0, 0}
	else 
		return {0, 0, 0, 0, rhoStar, 0, 0, 0}
	end
end

function Problem.Solver.MomEq.BC:WallsV(pos, initPos, states, t)
	return {0, 0, 0}
end
