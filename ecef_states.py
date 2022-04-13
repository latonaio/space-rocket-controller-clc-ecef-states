# TO-DO: Change microservice name
# TO-DO: Containerize to microservice
import numpy as np

# This microservice converts ECI coordinates to longitude, latitude, and altitude
# Requires: Initialization, instance of metadata distributor
# Required by: Environment model, direction cosines
# Outputs: longitude, latitude, altitude
class ECEFState():
  def __init__(self, metadata_distributor):
    self.metadata_distributor = metadata_distributor
    self.constants = self.metadata_distributor.constants
    self.xyz_s = self.metadata_distributor.get_var('xyz_s')
  
  def eci_to_lla(self, t):
    # position is instaneous (shape === (3,1))
    # t is time, to be decided
    x_s,y_s,z_s = self.xyz_s
    theta = np.arctan(self.constants.A/self.constants.B*z_s/(x_s**2+y_s**2)**0.5)
    # Longitude
    lamb = np.arctan(y_s - x_s, self.constants.OMEGA*t)
    # Latitude
    eta = np.arctan(((1-self.constants.E**2)*z_s+self.constants.E**2*self.constants.B*np.sin(theta)**3)/(
                    (1-self.constants.E**2)*((x_s**2+y_s**2)**0.5-self.constants.E**2*self.constants.A*np.cos(theta)**3)))
    # Altitude
    h = (x_s**2+y_s**2)**0.5/np.cos(eta) - self.constants.A/(1-self.constants.E**2*np.sin(eta)**2)
    self.metadata_distributor.set({'lamb': lamb, 'eta': eta, 'h': h})
    return lamb, eta, h
