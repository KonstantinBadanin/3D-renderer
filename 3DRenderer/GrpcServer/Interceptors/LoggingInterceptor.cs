using Grpc.Core;
using Grpc.Core.Interceptors;

namespace GrpcServer.Interceptors
{
    public class LoggingInterceptor: Interceptor
    {
        private ILogger<LoggingInterceptor> _logger;

        public LoggingInterceptor(ILogger<LoggingInterceptor> logger)
        {
            _logger = logger;
        }

        public override async Task DuplexStreamingServerHandler<TRequest, TResponse>(IAsyncStreamReader<TRequest> requestStream, IServerStreamWriter<TResponse> responseStream, ServerCallContext context, DuplexStreamingServerMethod<TRequest, TResponse> continuation)
        {
            _logger.Log(LogLevel.Information, $"Start - {context.Method}");

            await base.DuplexStreamingServerHandler(requestStream, responseStream, context, continuation);

            _logger.Log(LogLevel.Information, $"End - {context.Method}");

            await base.DuplexStreamingServerHandler(requestStream, responseStream, context, continuation);
        }
    }
}
