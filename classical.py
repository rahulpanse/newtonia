import matplotlib.pyplot as plt
import numpy as np
import math
import tkinter as tk
from tkinter import ttk, messagebox
import inspect
import textwrap

# Author: Rahul Panse
# Date: 10th September, 2025 
# Version: 1.0

"""Classical Mechanics Module
Provides functions to compute and plot displacement, velocity, and acceleration over time.
Other classical graphs can be added in the future.

TODO LIST:
- Add projectile motion calculations and plots.
- Add circular motion calculations and plots.
- Add harmonic motion calculations and plots.
-Add more features as needed."""

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
    def trig_function():
        function = input("Which function do you want to use? (sin, cos, tan): ").strip().lower()
        angle = float(input("Enter the angle in degrees: "))
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

# class DifferentialEquations:
#     def velocity():
#         highest_power = int(input("Enter the highest power of the differential: "))
#         no_of_terms = highest_power + 1
# 

def gui_help_menu():
    """
    Elegant 'website-like' Tkinter GUI Help Menu with collapsible sidebar.
    - Accordion-style navigation for Mechanics & Thermodynamics
    - Click a function -> opens modal with full documentation
    """

    root = tk.Tk()
    root.title("📘 Classical Mechanics & Thermodynamics Help")
    root.geometry("1000x650")
    root.configure(bg="#1e1e1e")

    # ----- STYLES -----
    MAIN_BG = "#1e1e1e"
    SIDEBAR_BG = "#252526"
    CONTENT_BG = "#2d2d30"
    HIGHLIGHT = "#00ff99"
    TEXT_COLOR = "#dddddd"
    FONT = ("Segoe UI", 11)

    # ----- SIDEBAR -----
    sidebar = tk.Frame(root, bg=SIDEBAR_BG, width=280)
    sidebar.pack(side="left", fill="y")

    title_label = tk.Label(
        sidebar,
        text="📘 Module Help",
        font=("Segoe UI", 16, "bold"),
        bg=SIDEBAR_BG,
        fg=HIGHLIGHT,
        pady=20
    )
    title_label.pack()

    # ----- MAIN CONTENT -----
    content_frame = tk.Frame(root, bg=CONTENT_BG)
    content_frame.pack(side="right", fill="both", expand=True)

    content_title = tk.Label(
        content_frame,
        text="Welcome!",
        font=("Segoe UI", 18, "bold"),
        bg=CONTENT_BG,
        fg=HIGHLIGHT
    )
    content_title.pack(pady=20)

    content_desc = tk.Label(
        content_frame,
        text="Expand a section on the left and click a function to view documentation.",
        font=FONT,
        bg=CONTENT_BG,
        fg=TEXT_COLOR,
        wraplength=700,
        justify="left"
    )
    content_desc.pack(pady=10)

    # ----- MODAL FUNCTION -----
    def show_modal(func_name, func_obj):
        doc = inspect.getdoc(func_obj)
        if not doc:
            doc = "No documentation available."

        modal = tk.Toplevel(root)
        modal.title(f"{func_name}() Documentation")
        modal.geometry("600x400")
        modal.configure(bg=CONTENT_BG)

        modal_title = tk.Label(
            modal,
            text=f"{func_name}()",
            font=("Segoe UI", 14, "bold"),
            bg=CONTENT_BG,
            fg=HIGHLIGHT
        )
        modal_title.pack(pady=10)

        doc_box = tk.Text(
            modal,
            wrap="word",
            bg=MAIN_BG,
            fg=TEXT_COLOR,
            insertbackground="white",
            font=("Consolas", 11),
            relief="flat"
        )
        doc_box.insert("1.0", textwrap.dedent(doc))
        doc_box.config(state="disabled")
        doc_box.pack(fill="both", expand=True, padx=15, pady=15)

        close_btn = tk.Button(
            modal,
            text="Close",
            command=modal.destroy,
            bg=HIGHLIGHT,
            fg="black",
            font=("Segoe UI", 10, "bold"),
            relief="flat",
            padx=10, pady=5
        )
        close_btn.pack(pady=10)

        modal.transient(root)
        modal.grab_set()
        root.wait_window(modal)

    # ----- MODULES -----
    modules = {
        "Mechanics ⚙️": Mechanics,
        "Thermodynamics 🌡️": Thermodynamics
    }

    # ----- COLLAPSIBLE SIDEBAR -----
    def toggle_section(frame, button):
        if frame.winfo_viewable():
            frame.pack_forget()
            button.config(text=button.cget("text").replace("▼", "▶"))
        else:
            frame.pack(fill="x", padx=20, pady=5)
            button.config(text=button.cget("text").replace("▶", "▼"))

    for module_name, module_class in modules.items():
        # Expand/Collapse Button
        section_btn = tk.Button(
            sidebar,
            text=f"▶ {module_name}",
            font=("Segoe UI", 13, "bold"),
            bg=SIDEBAR_BG,
            fg="#00aaff",
            activebackground="#333",
            activeforeground=HIGHLIGHT,
            relief="flat",
            anchor="w"
        )
        section_btn.pack(fill="x", padx=15, pady=5)

        # Frame that holds function buttons
        func_frame = tk.Frame(sidebar, bg=SIDEBAR_BG)

        for func_name, func_obj in inspect.getmembers(module_class, inspect.isfunction):
            btn = tk.Button(
                func_frame,
                text=f"{func_name}()",
                font=FONT,
                bg=SIDEBAR_BG,
                fg=TEXT_COLOR,
                activebackground="#333",
                activeforeground=HIGHLIGHT,
                relief="flat",
                anchor="w",
                command=lambda n=func_name, o=func_obj: show_modal(n, o)
            )
            btn.pack(fill="x", padx=25, pady=2)

        # Bind toggle behavior
        section_btn.config(command=lambda f=func_frame, b=section_btn: toggle_section(f, b))

    root.mainloop()



    
if __name__ == "__main__":
    gui_help_menu()