// See https://aka.ms/new-console-template for more information
using Grpc.Core;
using Grpc.Net.Client;
using GrpcServer;

using (var channel = GrpcChannel.ForAddress("https://localhost:7130"))
{
    var client = new Greeter.GreeterClient(channel);

    using (var call = client.SayHelloStream())
    {
        var readTask = Task.Run(async () =>
        {
            await foreach (var response in call.ResponseStream.ReadAllAsync())
            {
                Console.WriteLine(response.Message);
            }
        });

        while (true)
        {
            var result = Console.ReadLine();
            if (string.IsNullOrWhiteSpace(result))
            {
                break;
            }

            await call.RequestStream.WriteAsync(new HelloRequest() { Name=result});
        }

        await call.RequestStream.CompleteAsync();
        await readTask;
    }
    //var reply = await client.SayHelloAsync(new HelloRequest() { Name = "World" });
    //Console.WriteLine("Hello " + reply.Message);

    Console.ReadKey();
}
return await new Task<int>(() => 0);