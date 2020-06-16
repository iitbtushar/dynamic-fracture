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
  file = "mesh.data";
};

// define the tolerance and max number of iteration
// for the AM solver (implemented in the ForkModule)
// AM_max_iteration = 1 => just one iteration as in Fenics code of Jian

amImpExp = 
{
  AM_tolerance     = 1e-5;
  AM_max_iteration = 500;
  convergence      = "Damage";
  Carlsson         = false;
};


chain1 = 
{
  // note that: for this problem, this is just an empty file
  // as we do not need Dirichlet BCs for the phase field
  // you cannot remove the input so we need just an empty file
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
         elements = "notch";
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
               young = 3.09e3;
               ft    = 75.0;
               Gc  = 0.3;
               l0  = 0.5;
            };
          };
          
          writeCracks=
          {
            file     = "$(CASE_NAME)";
            interval = 1;
            x0       = 10.0;
            y0       = 80.;
            phiMax   = 0.9;
          };

        };
      };
    };
  };
  
  solver = "Nonlin"
  {
      precision  = 1.0e-6;
      //lineSearch = true;
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
        
          density   = 1180e-12;
	  E = 3.09e3; // to compute CFL time
          damThreshold = 0.0;//damage limit as a switch betn static and dynamic solver
          totalDisp  = 0.20;
    
          material = ["mat1"]; 
          mat1 =  "IsoPhaseField"
          {
            elements = "bulk";
            dim    = 2;
            state  = "PLANE_STRESS";
            thickness = 1.0;
    
            E     = 3.09e3;
            nu    = 0.35;
            
            degradation=  "WuDegradation"
            {
               softening="linear";
               young = 3.09e3;
               ft    = 75.0;
               Gc  = 0.3;
               l0  = 0.5;
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
          increment  = 0.10;          
	  velocity   = 0.;
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
          //lineSearch = true;
         
          solver = "CG"
          {
              precision = 1.0e-7;
              precon.type="ILUd";
              useThreads=true;
          };
       };

      exp =  
      {
        deltaTime = 2e-8;
      };
    };

    vtk = "VTKOutput"
    {
        prefix   = "vtu/$(CASE_NAME)";
        format   = "binary"; // binary and appended 
        elements = "bulk";
        interval = 1000;
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



