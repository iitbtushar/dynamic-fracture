mpart.overlap=1;

log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  pause = 0.;
  runWhile = "i < 0";
};
  
input =
{
  file = "d30.data";
};

// define the tolerance and max number of iteration

amImpExp = 
{
  AM_tolerance     = 1e-5;
  AM_max_iteration = 500;
  convergence      = "Damage";
  Carlsson         = true;
};


chain1 = 
{
  input =
  {
    file = "phi.data";
  };

  model = "MP" // for parallel calculations
  {
    model = "Matrix"
    {
      matrix.type = "FEM";
      matrix.symmetric = true;

    model       = "Multi"
    {
      models = ["dam","fixDam"];

      fixDam = "FixDomainPhi"
      {
         elements = "hinge";
	 Phi0 = 0.0;
      };
    
      dam =
      {
        type     = "VariationalDamPhase";
          elements = "bulk";
          elemType = "Lagrange";
        
          shape =
          {
            type      = "Triangle3";
            intScheme = "Gauss1";
          };

	  thickness = 1.0;
    
          geom  = "Wu";
          material  = ["mat"];

          mat = 
          {
            elements = "bulk";
            degradation=  "WuDegradation"
            {
               softening="linear";
               young = 4.39e3;
               ft    = 40.0;
               Gc  = 0.0933;
               l0  = 0.1;
            };
          };
          
          writeCracks=
          {
            file     = "$(CASE_NAME)";
            interval = 1;
            x0       = 8.;
            y0       = 60.;
            phiMax   = 0.9;
          };

        };
      };
    };
  };
  
  solver = "Nonlin"
  {
      precision  = 1.0e-6;
      bounds     = ["b1"];

      solver = "CG"
      {
          precision = 1.0e-7;
          precon.type="ILUd";
          useThreads=true;
      };
      
      b1 = 
      {
         dofType    = "pf";
         lowerBound = 0.;
         upperBound = 1.;
      };
   };
  
  extraModules =
  {
    modules =["energies"];
    
    energies = "Sample"
    {
      file       = "$(CASE_NAME)-energies.dat";
      sampleWhen = "i % 2"; // every 4 steps
      header     = "Time | surface energy | strain energy | kinetic energy";
      dataSets   = [ "t","surfaceEnergy","strainEnergy", "kineticEnergy"];
    };
    

    view = "FemView"
    {
      window =
      {
        height = 300;
        width  = 900;
      };
    
      snapFile = "$(CASE_NAME)%2i.png";
      configFile  = "$(CASE_NAME).view";
    
      dataSets = [ "disp" ];
    
      disp =
      {
        type   = "Vector";
        vector = "state";
      };
    
      mesh =
      {
        //showFaces = false;
    
        plugins = "colors";
        colors =
        {
          type = "MeshColorView";
          data = "disp";
          palette   = "custom";
          autoScale = true;
        };
      };
    };
  };
};

chain2 = 
{
  input =
  {
    file = "u.data";
  };

  model = "MP"
  {
    model = "Matrix"
    {
      matrix.type = "FEM";
      matrix.symmetric = true;
    
      model       = "Multi"
      {
        models = [ "bulk", "dirichlet", "lodi" ];

        lodi =  "Lodi"
        {
          group = "top";
        };
    
        bulk = "VariationalDamMech"
        {
          elements = "bulk";
          elemType = "Lagrange";
    
          shape =
          {
            type      = "Triangle3";
            intScheme = "Gauss1";
          };
        
          density   = 1170e-12;
	  E = 4.39e3; // to compute CFL time
          damThreshold = 0.65;//damage limit as a switch betn static and dynamic solver
          totalDisp  = 0.0; //for ref
	  reduceMass = false;
	  largeDeform = false;
    
          material = ["mat1"]; 
          mat1 =  "IsoPhaseField"
          {
            elements = "bulk";
            dim    = 2;
            state  = "PLANE_STRESS";
            thickness = 1.0;
    
            E     = 4.39e3;
            nu    = 0.37;
            
            degradation=  "WuDegradation"
            {
               softening="linear";
               young = 4.39e3;
               ft    = 40.0;
               Gc  = 0.0933;
               l0  = 0.1;
            };
            
            drivingForce = 
            {
              type = "Rankine";
            };
          };
        };
        
        dirichlet = "ImpExpDispControl" 
        {
          nodes1     = "top";
          nodes2     = "bottom";
          dof        = "dy";
          increment  = 1e-2;          
	  velocity   = 1e-2;
        };
      };
    };
  }; // end of MP model
  
  

  extraModules = 
  {
    modules = ["solver", "vtk", "sample"];

    solver = "ImpExpSolver"
    {
      nonLin = 
      {
          precision  = 1.0e-6;         
          solver = "CG"
          {
              precision = 1.0e-7;
              precon.type="ILUd";
              useThreads=true;
          };
       };

      exp =  
      {
        deltaTime = 3e-9;
	reduceMass = false;
      };
    };

    vtk = "VTKOutput"
    {
        prefix   = "vtu/$(CASE_NAME)";
        format   = "binary"; // binary and appended 
        elements = "bulk";
        interval = 5000;
        data     = "stress | damage";
    };

    sample = "Sample"
    {
      file       = "$(CASE_NAME).dat";
      header     = "   t |  uy[mm] | fy[N] ";
      dataSets   = ["t", " chain2.model.model.model.lodi.disp[1]"," chain2.model.model.model.lodi.load[1]"];    
    };

  };
};



