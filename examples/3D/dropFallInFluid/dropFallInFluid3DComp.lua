Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 1,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.05,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.35,
		boundingBox = {-1, -1, -1, 1, 1, 100},
		addOnFS = true,
		deleteFlyingNodes = false,
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/3D/dropFallInFluid/geometry.msh"
	},
	
	Extractors = {
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
		gamma = 0,
		K0 = 220000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		WallsFixed = true,
	},
	
	Solver = {
	    id = "CDS_dpdt",
		adaptDT = true,
		maxDT = 0.002,
		initialDT = 1e-8,
		securityCoeff = 0.5,
		
		MomEq = {
			bodyForce = {0, 0, -9.81},
			BC = {
			
			}
		},
		
		ContEq = {
			stabilization = "Meduri",
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
	local z0 = 0.3
	
	if(pos[3] <= z0) then
	    local rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[3]) + 1)^(1/(K0p - 1))
		local p = K0/K0p*((rho/rhoStar)^K0p - 1)
		return {0, 0, 0, p, rho, 0, 0, 0}
	else 
		return {0, 0, 0, 0, rhoStar, 0, 0, 0}
	end
end

function Problem.Solver.MomEq.BC:WallsV(pos, t)
	return {0, 0, 0}
end
