import unittest,math
from continuation_checks import check_endpoints,check_pulse
class RestartChecks(unittest.TestCase):
 def test_existing_checkpoint_schedule(self):
  self.assertEqual(check_endpoints(list(range(0,12537,16)),783,12536,16),783)
 def test_resume_schedule_starts_after_saved_step(self):
  self.assertEqual(check_endpoints([12544],1,12544,16,start=12536),1)
 def test_resume_missing_teacher_rejected(self):
  with self.assertRaises(AssertionError):check_endpoints([],1,12544,16,start=12536)
 def test_resume_extra_teacher_rejected(self):
  with self.assertRaises(AssertionError):check_endpoints([12528,12544],1,12544,16,start=12536)
 def test_resume_wrong_svd_count_rejected(self):
  with self.assertRaises(AssertionError):check_endpoints([12544],2,12544,16,start=12536)
 def test_resume_uses_absolute_pulse_time(self):
  dt=.03;omega=.0569614493722652368;tw=882.447392251770225
  def ac(t):
   tt=t-tw/2
   return -5.338e-9*1e6/omega*math.cos(math.pi*tt/tw)**2*math.sin(omega*tt)
  rows=[]
  for i in range(12537,12545):
   a=ac(i*dt);e=-(a-ac((i-1)*dt))/dt
   rows.append([i*dt,0.,0.,a,0.,0.,e,0.,0.,a,0.,0.,e,0.,0.,0.])
  self.assertLess(check_pulse(rows,start=12536),1e-14)
  rows[0][3]=0.
  with self.assertRaises(AssertionError):check_pulse(rows,start=12536)
if __name__=='__main__':unittest.main()
