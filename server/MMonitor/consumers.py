from channels.generic.websocket import AsyncWebsocketConsumer
import json

class MyConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        await self.accept()

    async def disconnect(self, close_code):
        pass

    async def receive(self, text_data):
        text_data_json = json.loads(text_data)
        message = text_data_json['message']

        # Do something with the message
        # Perhaps send it to other connected clients, save to DB, etc.

        await self.send(text_data=json.dumps({ 
            'message': message
        }))
