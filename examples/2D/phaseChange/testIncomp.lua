Problem = {
    id = "Boussinesq",
	simulationTime = 1.5,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.0001,
		alpha = 1.2,
		omega = 0.85,
		gamma = 0.6,
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
			whatToWrite = {"T", "ke", "p", "velocity"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1.3e-3,
		rho = 2385,
		k = 94.03,
		alpha = 1.17e-4,
		Tr = 660 + 273.15,
		cv = 1080,
		gamma = 0.812,
		DgammaDT = -2e-4,
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
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.00001,
		initialDT = 0.00001,
		solveHeatFirst = true,
		maxRemeshDT = -1,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, -9.81},
			residual = "Ax_f",
			BC = {

			}
		},
		
		HeatEq = {
			minRes = 1e-6,
			maxIter = 10,
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0, 300}
end

--[[function Problem.Solver.HeatEq.BC:TopQ(pos, initPos, t) 
	return {0}
end]]--

function Problem.Solver.MomContEq.BC:TopV(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:BottomQ(pos, t) 
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:BottomV(pos, t) 
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

function Problem.Solver.MomContEq.BC:LeftV(pos, t) 
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

function Problem.Solver.MomContEq.BC:RightV(pos, t) 
	return {0, 0}
end
