import matplotlib.pyplot as plt
import numpy as np
import math
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
import inspect
import textwrap
import ttkbootstrap as tb
import sys

# Author: Rahul Panse
# Date: 10th September, 2025 
# Version: 1.0

"""Classical Mechanics Module
Provides functions to compute and plot displacement, velocity, and acceleration over time.
Other classical graphs can be added in the future.

TODO LIST:
- Add harmonic motion calculations and plots.
- Add more features as needed."""

class Mechanics:
    def disp_time(initial_displacement, initial_velocity, acceleration, time):
        """
        Calculate displacement over time given initial conditions and constant acceleration.
        
        Parameters:
        initial_displacement (float): Initial displacement in meters.
        initial_velocity (float): Initial velocity in meters per second.
        acceleration (float): Constant acceleration in meters per second squared.
        time (numpy array): Array of time values in seconds.
        
        Returns:
        numpy array: Displacement values corresponding to the time array.
        """
        time = np.linspace(0, time, 100)  # seconds

        # Displacement calculation: s = s0 + v0*t + 0.5*a*t^2
        displacement = initial_displacement + initial_velocity * time + 0.5 * acceleration * time**2

        # Plotting
        plt.figure(figsize=(8, 5))
        plt.plot(time, displacement, label='Displacement')
        plt.xlabel('Time (s)')
        plt.ylabel('Displacement (m)')
        plt.title('Displacement-Time Graph')
        plt.gcf().canvas.manager.set_window_title('Displacement-Time Computation')  # Set window title
        plt.legend()
        plt.grid(True)
        plt.show()
    def velocity_time(initial_velocity, acceleration, time):
        """
        Calculate velocity over time given initial velocity and constant acceleration.
        
        Parameters:
        initial_velocity (float): Initial velocity in meters per second.
        acceleration (float): Constant acceleration in meters per second squared.
        time (numpy array): Array of time values in seconds.
        
        Returns:
        numpy array: Velocity values corresponding to the time array.
        """
        time = np.linspace(0, time, 100)  # seconds

        # Velocity calculation: v = v0 + a*t
        velocity = initial_velocity + acceleration * time

        # Plotting
        plt.figure(figsize=(8, 5))
        plt.plot(time, velocity, label='Velocity', color='orange')
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity (m/s)')
        plt.title('Velocity-Time Graph')
        plt.gcf().canvas.manager.set_window_title('Velocity-Time Computation')  # Set window title
        plt.legend()
        plt.grid(True)
        plt.show()
    def acceleration_time(initial_velocity, final_velocity, time):
        """
        Calculate acceleration over time given initial and final velocities.
        
        Parameters:
        initial_velocity (float): Initial velocity in meters per second.
        final_velocity (float): Final velocity in meters per second.
        time (numpy array): Array of time values in seconds.
        
        Returns:
        numpy array: Acceleration values corresponding to the time array.
        
        time = np.linspace(1,time, time + 10)
        acceleration = np.full_like(time, (final_velocity - initial_velocity) / time)
        plt.figure(figsize=(8, 5))
        plt.gcf().canvas.manager.set_window_title('Acceleration-Time Simulation')  # Set window
        plt.xlabel('Time (s)')
        plt.ylabel('Acceleration (m/s²)')
        plt.title('Acceleration-Time Graph')
        plt.legend()
        #plt.grid(True)
        plt.plot(time, acceleration, label='Acceleration', color='green')
        plt.show()""" 
        pass
    def force_displacement(force, displacement):
        """
        Calculate work done given force and displacement.
        
        Parameters:
        force (float): Force in Newtons.
        displacement (float): Displacement in meters.
        
        Returns:
        float: Work done in Joules.
        """
        work_done = force * displacement  # Work done calculation: W = F * d
        print(f"Work done: {work_done} Joules")
    def gravitational_force(mass1, mass2, distance):
        """
        Calculate gravitational force between two masses.
        
        Parameters:
        mass1 (float): Mass of the first object in kilograms.
        mass2 (float): Mass of the second object in kilograms.
        distance (float): Distance between the centers of the two masses in meters.
        
        Returns:
        float: Gravitational force in Newtons.
        """
        G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
        force = G * (mass1 * mass2) / distance**2  # Gravitational force calculation: F = G * (m1*m2)/r^2
        print(f"Gravitational Force: {force} Newtons")
    def compute_angle():
        numerical_value = float(input("Enter the numerical value: "))
        get_func = input("Which function do you want to use? (sin, cos, tan): ").strip().lower()
        if get_func == "sin":
            angle = math.sinh(numerical_value)
            print(f"Your angle is: {math.degrees(angle)} degrees")
        elif get_func == "cos":
            angle = math.cosh(numerical_value)
            print(f"Your angle is: {math.degrees(angle)} degrees")
        elif get_func == "tan":
            angle = math.tanh(numerical_value)
            print(f"Your angle is: {math.degrees(angle)} degrees")
        else:
            print("Invalid function. Please choose 'sin', 'cos', or 'tan'.")
            return
    def trig_function(trig_func, angle):
        function = trig_func.strip().lower()
        angle_rad = math.radians(angle)  # Convert angle to radians
        if function == "sin":
            result = math.sin(angle_rad)
            print(f"sin({angle} degrees) = {result}")
        elif function == "cos": 
            result = math.cos(angle_rad)
            print(f"cos({angle} degrees) = {result}")   
        elif function == "tan":
            result = math.tan(angle_rad)
            print(f"tan({angle} degrees) = {result}")
        else:
            print("Invalid function. Please choose 'sin', 'cos', or 'tan'.")
            return
    class ProjectileMotion:
        def maximum_height(initial_velocity, angle):
            """
            Calculate the maximum height of a projectile.
            
            Parameters:
            initial_velocity (float): Initial velocity in meters per second.
            angle (float): Launch angle in degrees.
            
            Returns:
            float: Maximum height in meters.
            """
            g = 9.81  # Acceleration due to gravity in m/s²
            angle_rad = math.radians(angle)  # Convert angle to radians
            max_height = (initial_velocity**2 * math.sin(angle_rad)**2) / (2 * g)
            print(f"Maximum Height: {max_height} meters")
            return max_height
        def range_of_projectile(initial_velocity, angle):
            """ 
            Calculate the range of a projectile.

            Parameters:
            initial_velocity (float): Initial velocity in meters per second.
            angle (float): Launch angle in degrees.
            
            Returns:
            float: Range in meters.
            """
            g = 9.81
            angle_rad = math.radians(angle)
            range_proj = (initial_velocity**2 * math.sin(2 * angle_rad)) / g
            print(f"Range of Projectile: {range_proj} meters")
            return range_proj
        def time_of_flight(initial_velocity, angle):
            """
            Calculate the time of flight of a projectile.
            
            Parameters:
            initial_velocity (float): Initial velocity in meters per second.
            angle (float): Launch angle in degrees.
            
            Returns:
            float: Time of flight in seconds.
            """
            g = 9.81
            angle_rad = math.radians(angle)
            time_flight = (2 * initial_velocity * math.sin(angle_rad)) / g
            print(f"Time of Flight: {time_flight} seconds")
            return time_flight
        def get_trajectory(initial_velocity, angle):
            """
            Calculate the trajectory of a projectile and plot it.
            
            Parameters:
            initial_velocity (float): Initial velocity in meters per second.
            angle (float): Launch angle in degrees.
            """
            g = 9.81
            angle_rad = math.radians(angle)
            time_of_flight = (2 * initial_velocity * math.sin(angle_rad)) / g
            t = np.linspace(0, time_of_flight, num=500)
            x = initial_velocity * np.cos(angle_rad) * t
            y = initial_velocity * np.sin(angle_rad) * t - 0.5 * g * t**2

            plt.figure(figsize=(10, 5))
            plt.plot(x, y)
            plt.title('Projectile Motion Trajectory')
            plt.xlabel('Distance (m)')
            plt.ylabel('Height (m)')
            plt.gcf().canvas.manager.set_window_title('Projectile Motion Trajectory')  # Set window title
            plt.ylim(bottom=0)
            plt.grid()
            plt.show()

class Thermodynamics:
    def ideal_gas_law(pressure=None, volume=None, moles=None, temperature=None):
        """
        Calculate the missing variable in the Ideal Gas Law equation: PV = nRT
        
        Parameters:
        pressure (float): Pressure in Pascals (Pa).
        volume (float): Volume in cubic meters (m^3).
        moles (float): Number of moles (n).
        temperature (float): Temperature in Kelvin (K).
        
        Returns:
        float: The calculated missing variable.
        """
        R = 8.314  # Ideal gas constant in J/(mol·K)
        
        if pressure is None:
            if volume is None or moles is None or temperature is None:
                raise ValueError("To calculate pressure, volume, moles, and temperature must be provided.")
            pressure = (moles * R * temperature) / volume
            print(f"Calculated Pressure: {pressure} Pa")
            return pressure
        
        elif volume is None:
            if pressure is None or moles is None or temperature is None:
                raise ValueError("To calculate volume, pressure, moles, and temperature must be provided.")
            volume = (moles * R * temperature) / pressure
            print(f"Calculated Volume: {volume} m^3")
            return volume
        
        elif moles is None:
            if pressure is None or volume is None or temperature is None:
                raise ValueError("To calculate moles, pressure, volume, and temperature must be provided.")
            moles = (pressure * volume) / (R * temperature)
            print(f"Calculated Moles: {moles} mol")
            return moles
        
        elif temperature is None:
            if pressure is None or volume is None or moles is None:
                raise ValueError("To calculate temperature, pressure, volume, and moles must be provided.")
            temperature = (pressure * volume) / (moles * R)
            print(f"Calculated Temperature: {temperature} K")
            return temperature
        
        else:
            raise ValueError("One of the parameters must be None to calculate it.")
    def flot(delta_u, q, w):
        """Calculates the First Law of Thermodynamics: ΔU = Q - W
        where ΔU is the change in internal energy, Q is the heat added to the system, and W is the work done by the system.
        Parameters are to be taken as user input.
        """
        
        # Calculating the missing variable
        if delta_u == None:
            if q is None or w is None:
                raise ValueError("To calculate ΔU, Q and W must be provided.")
            change_in_internal_energy = q - w
            print(f"Calculated Change in Internal Energy (ΔU): {change_in_internal_energy}")
        elif q == None:
            if delta_u is None or w is None:
                raise ValueError("To calculate Q, ΔU and W must be provided.")
            heat_added = delta_u + w
            print(f"Calculated Heat Added (Q): {heat_added}")
        elif w == None:
            if delta_u is None or q is None:
                raise ValueError("To calculate W, ΔU and Q must be provided.")
            work_done = q - delta_u
            print(f"Calculated Work Done (W): {work_done}")
        else:
            raise ValueError("[!] One of the parameters mus be None to calculate it.")
    def slot():
        """Calculates the second law of thermodynamics using the formula: ΔS = Q/T
        where ΔS is the change in entropy, Q is the heat added to the system, and T is the absolute temperature in Kelvin.
        Parameters are to be taken as user input.
        
        There are many types of entropy calculations, I will be adding them in the quantum mechanics module.
        von Neumann entropy, Shannon entropy, etc. This is just a basic one. I will also be using S = kb * ln(omega) in this particular module. 
        Further, the second formula will always contain the Boltzmann constant, kb = 1.38e-23 J/K and omega, which is the number of microstates.
        """
        heat_added = float(input("Enter the heat added to the system (Q) in Joules: "))
        temperature = float(input("Enter the absolute temperature (T) in Kelvin: "))
        
        if temperature <= 0:
            raise ValueError("Temperature must be greater than 0 Kelvin.")
        
        change_in_entropy = heat_added / temperature
        print(f"Calculated Change in Entropy (ΔS): {change_in_entropy} J/K")
        """Another formula for entropy is S = kb * ln(omega)
        where kb is the Boltzmann constant and omega is the number of microstates. Section 2"""
        boltzmann_constant = 1.38e-23  # J/K
        microstates = float(input("Enter the number of microstates (omega): "))
        if microstates <= 0:
            print(f"[!] The number of microstates must be greater than 0. \nError: {microstates} is not valid.\nException: {ValueError("Invalid number of microstates.")}")
        else:
            entropy = boltzmann_constant * math.log(microstates)
            print("*******************Formula 2*******************")
            print(f"Calculated Entropy (S) using S = kb * ln(omega): {entropy} J/K")

class Rotational:
    def angular_velocity(radius, tangential_velocity):
        """
        Calculate the angular velocity of a rotating object.
        
        Parameters:
        radius (float): Radius in meters.
        tangential_velocity (float): Tangential velocity in meters per second.
        
        Returns:
        float: Angular velocity in radians per second.
        """
        angular_vel = tangential_velocity / radius
        print(f"Angular Velocity: {angular_vel} radians/second")
        return angular_vel
    def angular_momentum(moment_of_inertia, angular_velocity):
        """ n
        Calculate the angular momentum of a rotating object.
        
        Parameters:
        moment_of_inertia (float): Moment of inertia in kg·m².
        angular_velocity (float): Angular velocity in radians per second.
        
        Returns:
        float: Angular momentum in kg·m²/s.
        """
        angular_mom = moment_of_inertia * angular_velocity
        print(f"Angular Momentum: {angular_mom} kg·m²/s")
        return angular_mom
    def rotational_kinetic_energy(moment_of_inertia, angular_velocity):
        """
        Calculate the rotational kinetic energy of a rotating object.
        
        Parameters:
        moment_of_inertia (float): Moment of inertia in kg·m².
        angular_velocity (float): Angular velocity in radians per second.
        
        Returns:
        float: Rotational kinetic energy in Joules.
        """
        rotational_ke = 0.5 * moment_of_inertia * angular_velocity**2
        print(f"Rotational Kinetic Energy: {rotational_ke} Joules")
        return rotational_ke
    def torque(force, radius, angle):
        """"The function calculates the torque using the formula
        τ = r * F * sin(θ)
        where τ is the torque, r is the radius (distance from the pivot point to the point where the force is applied)
        and theta is the angle between the force vector and the radius vector. Originally, the formula is 
        the cross-product of r and F.
        
        Params: 
        force (float): Force in Newtons.
        radius (float): Radius in meters.
        angle (float): Angle in degrees between the force vector and the radius vector.

        The cross product of any two vectors is given as A x B = |A|.|B|sin(theta)
        """
        torque = radius * force * math.sin(math.radians(angle))
        print(f"Torque: {torque} Newton-meters")
        return torque
    def mominertia_vs_mass():
        """This function plots the moment of inertia vs mass for a solid sphere and a solid cylinder."""
        mass = np.linspace(1, 100, 100)  # Mass in kg
        radius = 1  # Radius in meters (constant for this example)

        # Moment of inertia calculations
        moment_of_inertia_sphere = (2/5) * mass * radius**2
        moment_of_inertia_cylinder = (1/2) * mass * radius**2

        # Plotting
        plt.figure(figsize=(8, 5))
        plt.plot(mass, moment_of_inertia_sphere, label='Solid Sphere', color='blue')
        plt.plot(mass, moment_of_inertia_cylinder, label='Solid Cylinder', color='red')
        plt.xlabel('Mass (kg)')
        plt.ylabel('Moment of Inertia (kg·m²)')
        plt.title('Moment of Inertia vs Mass')
        plt.gcf().canvas.manager.set_window_title('Moment of Inertia vs Mass')  # Set window title
        plt.legend()
        plt.grid(True)
        plt.show()

class GUI:
    def help_gui():
        """
        Opens an interactive GUI documenting the contents of this module:
        - Collapsible class tabs for each top-level class
        - Collapsible function panels for each function inside a class
        - Shows function signature and docstring
        - Allows calling functions interactively (simple input parsing)

        Optimized with scrollable right details area so content never gets cut off.
        """

        # --- styling / theme (use ttkbootstrap if available) ---
        
        
        root = tb.Window(themename="solar")  
        style = tb.Style()
        THEMES = style.theme_names()
        using_ttkb = True


        root.title("Classical Module — Interactive Documentation")
        root.geometry("1100x750")  # slightly larger default

        # top frame
        top = ttk.Frame(root, padding=(12, 8, 12, 0))
        top.pack(fill="x")
        ttk.Label(top, text="Classical Module — Interactive Help", font=("Helvetica", 16, "bold")).pack(side="left")
        search_var = tk.StringVar()
        search_entry = ttk.Entry(top, width=36, textvariable=search_var)
        search_entry.pack(side="right", padx=(6,2))
        ttk.Label(top, text="Search:").pack(side="right")

        # main layout
        paned = ttk.Panedwindow(root, orient="horizontal")
        paned.pack(fill="both", expand=True, padx=12, pady=12)
        left_frame = ttk.Frame(paned, width=300)
        right_frame = ttk.Frame(paned)
        paned.add(left_frame, weight=1)
        paned.add(right_frame, weight=3)

        # --- LEFT SCROLL NAV ---
        left_canvas = tk.Canvas(left_frame, borderwidth=0, highlightthickness=0)
        left_scroll = ttk.Scrollbar(left_frame, orient="vertical", command=left_canvas.yview)
        left_inner = ttk.Frame(left_canvas)
        left_inner.bind("<Configure>", lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")))
        left_canvas.create_window((0,0), window=left_inner, anchor="nw")
        left_canvas.configure(yscrollcommand=left_scroll.set)
        left_canvas.pack(side="left", fill="both", expand=True)
        left_scroll.pack(side="right", fill="y")

        # --- RIGHT SCROLL DETAILS ---
        right_canvas = tk.Canvas(right_frame, borderwidth=0, highlightthickness=0)
        right_scroll = ttk.Scrollbar(right_frame, orient="vertical", command=right_canvas.yview)
        right_inner = ttk.Frame(right_canvas)
        right_inner.bind("<Configure>", lambda e: right_canvas.configure(scrollregion=right_canvas.bbox("all")))
        right_canvas.create_window((0,0), window=right_inner, anchor="nw")
        right_canvas.configure(yscrollcommand=right_scroll.set)
        right_canvas.pack(side="left", fill="both", expand=True)
        right_scroll.pack(side="right", fill="y")

        details_title = ttk.Label(right_inner, text="Select a function to view details", font=("Helvetica", 14, "bold"))
        details_title.pack(anchor="w", pady=(0,6))

        details_area = ttk.Frame(right_inner)
        details_area.pack(fill="both", expand=True)

        # --- Collapsible helper ---
        class Collapsible(ttk.Frame):
            def __init__(self, parent, text=""):
                super().__init__(parent)
                self._open = tk.BooleanVar(value=False)
                self.btn = ttk.Checkbutton(self, text=text, style="Toolbutton", variable=self._open, command=self._toggle)
                self.btn.pack(fill="x", padx=4, pady=2)
                self.container = ttk.Frame(self)
            def _toggle(self):
                if self._open.get():
                    self.container.pack(fill="both", expand=True, padx=8, pady=(0,6))
                else:
                    self.container.forget()

        # format helper
        def format_doc(obj):
            try: sig = str(inspect.signature(obj))
            except Exception: sig = "(signature unavailable)"
            doc = inspect.getdoc(obj) or "(no docstring)"
            return sig, textwrap.dedent(doc)

        # --- build function invoker ---
        def build_invoker(parent, func):
            frm = ttk.Frame(parent, padding=8)
            frm.pack(fill="x")
            sig, doc = format_doc(func)
            ttk.Label(frm, text=f"Signature: {sig}", font=("TkDefaultFont", 10, "italic")).pack(anchor="w")
            ttk.Separator(frm, orient="horizontal").pack(fill="x", pady=6)

            doc_text = tk.Text(frm, height=8, wrap="word")
            doc_text.insert("1.0", doc)
            doc_text.configure(state="disabled")
            doc_text.pack(fill="both", expand=False)

            # args
            params=[]
            try:
                for name, par in inspect.signature(func).parameters.items():
                    ttk.Label(frm, text=f"{name} =").pack(anchor="w", pady=(6,0))
                    ent=ttk.Entry(frm); ent.pack(fill="x"); params.append((name, ent))
            except Exception: pass

            outbox = tk.Text(frm, height=6, wrap="word"); outbox.pack(fill="both", pady=(8,0))

            def try_parse(s):
                if not s: return None
                s=s.strip()
                if s.lower()=="none": return None
                if s.lower() in ("true","false"): return s.lower()=="true"
                if s.startswith("[") and s.endswith("]"):
                    return [try_parse(x) for x in s[1:-1].split(",")]
                try:
                    return float(s) if "." in s else int(s)
                except: pass
                if s.startswith("np."):
                    try: return eval(s, {"np":np, "math":math})
                    except: pass
                return s

            def invoke():
                outbox.configure(state="normal"); outbox.delete("1.0","end")
                kwargs={}
                for n,e in params:
                    try: kwargs[n]=try_parse(e.get())
                    except Exception as ex: outbox.insert("end",f"Error parsing {n}: {ex}\n"); return
                try:
                    res=func(**kwargs) if kwargs else func()
                    outbox.insert("end", f"Result: {repr(res)}\n" if res is not None else "[Done] Function returned None\n")
                except Exception as ex:
                    outbox.insert("end", f"Exception: {ex}\n")
                outbox.configure(state="disabled")

            ttk.Button(frm,text="Run",command=invoke).pack(pady=6)
            return frm

        # --- introspect this module ---
        this_module = sys.modules[__name__]
        classes=[(n,o) for n,o in vars(this_module).items() if inspect.isclass(o) and o.__module__==this_module.__name__]

        nav_buttons=[]
        def show_details_for(func_obj, fullname):
            for w in details_area.winfo_children(): w.destroy()
            ttk.Label(details_area, text=fullname, font=("Helvetica", 13, "bold")).pack(anchor="w")
            build_invoker(details_area, func_obj)

        for cname,cls in sorted(classes,key=lambda x:x[0].lower()):
            col=Collapsible(left_inner,text=f"Class: {cname}")
            col.pack(fill="x", pady=4)
            cls_doc=(inspect.getdoc(cls) or "").split("\n")[0]
            if cls_doc: ttk.Label(col.container,text=cls_doc,wraplength=260,foreground="gray").pack(anchor="w", padx=6)
            members=[]
            for mname,mobj in sorted(vars(cls).items()):
                if (inspect.isfunction(mobj) or inspect.ismethod(mobj)) and not mname.startswith("_"):
                    members.append((mname,mobj))
                if inspect.isclass(mobj) and mobj.__module__==this_module.__name__:
                    for nm,no in sorted(vars(mobj).items()):
                        if (inspect.isfunction(no) or inspect.ismethod(no)) and not nm.startswith("_"):
                            members.append((f"{mname}.{nm}", no))
            for fname,fobj in members:
                fbtn=ttk.Button(col.container,text=fname,style="Toolbutton",
                    command=lambda fo=fobj,fn=f"{cname}.{fname}": show_details_for(fo,fn))
                fbtn.pack(fill="x", padx=6, pady=2)
                nav_buttons.append((f"{cname}.{fname}",fbtn))

        def do_search(*_):
            q=search_var.get().lower()
            for fullname,btn in nav_buttons:
                btn.pack_forget()
                if not q or q in fullname.lower():
                    btn.pack(fill="x", padx=6, pady=2)
        search_var.trace_add("write",do_search)

        root.protocol("WM_DELETE_WINDOW", lambda: root.destroy())
        root.mainloop()


if __name__ == "__main__":
    GUI.help_gui()