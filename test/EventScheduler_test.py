import unittest
import weakref
import _gfrd as mod
from gfrdbase import DomainEvent


#class Delegate(object):

#    def __init__(self, obj, method):
#        self.obj = weakref.proxy(obj)
#        self.method = method

#    def __call__(self, arg):
#        return self.method(self.obj, arg)



class EventSchedulerTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        scheduler = mod.EventScheduler()
        self.assertFalse(scheduler == None)

    def test_empty_state(self):
        scheduler = mod.EventScheduler()
        self.assertFalse(scheduler.size != 0)
        self.assertFalse(scheduler.time != 0.0)
        # what if getTopEvent() are called here?

    def test_one_event(self):
        scheduler = mod.EventScheduler()

        event = mod.PythonEvent(1.0, "TestEvent")
        id = scheduler.add(event)

        self.assertEqual(scheduler.time , 0.0)
        self.assertEqual(scheduler.top[1].time ,  1.0)
        self.assertEqual(scheduler.top[1].data , "TestEvent")
        self.assertEqual(scheduler.top[0] , id)

        self.assertEqual((id, event), scheduler.pop())
        self.assertEqual(scheduler.size , 0)
        self.assertEqual(scheduler.time , 1.0)


    def test_peek_second_event(self):

        scheduler = mod.EventScheduler()

        event1 = mod.PythonEvent(1.0, "Ev1")
        event2 = mod.PythonEvent(0.5, "Ev2")

        id1 = scheduler.add(event1)
        id2 = scheduler.add(event2)

        self.assertEqual(2, scheduler.size)

        second = scheduler.second

        self.assertEqual(1.0, second[1].time)
        self.assertEqual(id1, second[0])




if __name__ == "__main__":
    unittest.main()
