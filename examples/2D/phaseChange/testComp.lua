Problem = {
    id = "BoussinesqWC",
	simulationTime = 1.5,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.001,
		alpha = 1.2,
		omega = 0.65,
		gamma = 0.5,
		boundingBox = {-0.03, -0.002, 0.03, 100},
		deleteFlyingNodes = false,
		addOnFS = true,
		laplacianSmoothingBoundaries = false,
		exclusionZones = {},
		mshFile = "examples/2D/thermalConv/geometry2.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.002,
			whatToWrite = {"T", "ke", "p", "rho", "velocity"},
			writeAs = "NodesElements" 
		}
	},
	
	MeshSmoother = {
		--[[{
			kind = "GETMe",
			epsTol = 0.01,
			maxIter = 20
		}--]]
	},
	
	Material = {
		mu = 1.3e-3,
		k = 94.03,
		alpha = 1.17e-4,
		Tr = 660 + 273.15,
		cv = 1080,
		gamma = 0.812,
		DgammaDT = -2e-4,
		K0 = 22000000,
		K0p = 7.6,
		rhoStar = 2385,
		h = 0,
		epsRad = 0,
		Tinf = 300,
		Tm = 660 + 273.15,
		DT = 20,
		C = 1e6,
		eps = 1e-3,
		Lm = 396000
	},
	
	IC = {
		BottomFixed = true,
		LeftFixed = true,
		RightFixed = true
	},
	
	Solver = {
	    id = "CDS_dpdt",
		adaptDT = true,
		maxDT = 0.00001,
		initialDT = 1e-8,
		securityCoeff = 0.1,
		
		MomEq = {
			bodyForce = {0, -9.81},
			BC = {
			
			}
		},
		
		ContEq = {
			stabilization = "Meduri",
			BC = {

			}
		},
		
		HeatEq = {
			BC = {
				
			}
		}
	}
}

function Problem.IC:initStates(pos)
	local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local rhoStar = Problem.Material.rhoStar
	local g = -Problem.Solver.MomEq.bodyForce[2]
	local Tr = Problem.Material.Tr
	local z0 = 0.005

	if(pos[2] < z0) then
		local rho, p
		rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[2]) + 1)^(1/(K0p - 1))
		p = K0/K0p*((rho/rhoStar)^K0p - 1)
		
		return {0, 0, p, rho, 0, 0, 300}
	else
		return {0, 0, 0, rhoStar, 0, 0, 300}
	end
end

function Problem.Solver.MomEq.BC:TopV(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:BottomQ(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:FreeSurfaceQ(pos, t)
	if(pos[1] > 0.01 + 0.003 or pos[1] < 0.01 - 0.003) then
		return {0, 0}
	else
		local Q = -175000000*(1 - math.abs(pos[1] - 0.01)/0.003)*(1 + math.abs(pos[1] - 0.01)/0.003)
		local dt = 0.05
		if(t < dt) then
			return {0, (-2*Q*(t/dt)^3 + 3*Q*(t/dt)^2)}
		elseif(t > dt and t < 4*dt) then
			return {0, Q}
		elseif(t > 4*dt and t < 4*dt) then
			return {0, (2*Q*((t - 4*dt)/dt)^3 - 3*Q*((t - 4*dt)/dt)^2) + Q}
		else
			return {0, 0}
		end
	end
end

function Problem.Solver.MomEq.BC:BottomV(pos, t) 
	return {0, 0}
end

-- function Problem.Solver.HeatEq.BC:LeftT(pos, t) 
	-- local Tr = Problem.Material.Tr
	-- local Th = Tr + 50
	-- local DT = 0.5
	-- if(t > DT) then
		-- return {Th}
	-- else
		-- return {(-2*(Th - Tr)*(t/DT)^3 + 3*(Th - Tr)*(t/DT)^2) + Tr}
	-- end
-- end

function Problem.Solver.MomEq.BC:LeftV(pos, t) 
	return {0, 0}
end

-- function Problem.Solver.HeatEq.BC:RightT(pos, t)
	-- local Tr = Problem.Material.Tr 
	-- local Tc = Tr - 50
	-- local DT = 0.5
	-- if(t > DT) then
		-- return {Tc}
	-- else
		-- return {(-2*(Tc - Tr)*(t/DT)^3 + 3*(Tc - Tr)*(t/DT)^2) + Tr}
	-- end
-- end

function Problem.Solver.MomEq.BC:RightV(pos, t) 
	return {0, 0}
end

